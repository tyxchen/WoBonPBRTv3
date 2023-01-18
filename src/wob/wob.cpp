//
// Created by Terry Chen on 2022-10-21.
//

#include "core/camera.h"
#include "core/film.h"
#include "core/imageio.h"
#include "core/paramset.h"
#include "core/progressreporter.h"
#include "core/scene.h"

#include "wob/wob.h"

#include <fstream>
#include <iomanip>

using namespace pbrt;
using namespace pbrt_ext;

inline Float G(const Point3f &x, const Point3f &y) {
    return Inv4Pi / (x - y).Length();
}

inline Float H(const Point3f &x, const Point3f &y, const Normal3f &n_y) {
    auto r = (y - x).Length();
    return -Inv4Pi * Dot(y - x, n_y) / (r * r * r);
}

WoBIntegrator::WoBIntegrator(int maxDepth, std::shared_ptr<const Camera> camera,
                             std::shared_ptr<Sampler> sampler,
                             const pbrt::Bounds2i &pixelBounds, pbrt::Float rrThreshold,
                             const std::string &boundaryCond,
                             const std::string &domainType,
                             const std::string &colourmapParam) :
    pbrt::SamplerIntegrator(std::move(camera), sampler, pixelBounds),
    maxDepth(5), rrThreshold(rrThreshold)
{
    if (boundaryCond == "dirichlet") {
        cond = DIRICHLET;
    } else if (boundaryCond == "neumann") {
        cond = NEUMANN;
    } else if (boundaryCond == "robin") {
        cond = ROBIN;
    } else {
        throw std::logic_error("Bad boundary condition given.");
    }

    if (domainType == "interior") {
        domain = INTERIOR;
    } else if (domainType == "exterior") {
        domain = EXTERIOR;
    } else {
        throw std::logic_error("Bad domain type given.");
    }

    if (colourmapParam == "jet") {
        colourmap = tinycolormap::ColormapType::Jet;
    } else if (colourmapParam == "heat") {
        colourmap = tinycolormap::ColormapType::Heat;
    } else if (colourmapParam == "turbo") {
        colourmap = tinycolormap::ColormapType::Turbo;
    } else if (colourmapParam == "viridis") {
        colourmap = tinycolormap::ColormapType::Viridis;
    } else {
        throw std::logic_error("Unknown colourmap given.");
    }
}

void WoBIntegrator::Render(const pbrt::Scene &scene)
{
    this->pbrt::SamplerIntegrator::Render(scene);

    // Save final image after rendering
    auto width = camera->film->croppedPixelBounds.Diagonal().x;

    size_t i = 0;
    std::ofstream out_file(camera->film->filename + ".txt");
    out_file << std::setprecision(3);

    for (auto pixel : camera->film->croppedPixelBounds) {
        auto &pixel_buf = camera->film->GetPixel(pixel);
        Float result[3];
        XYZToRGB(pixel_buf.xyz, result);
        auto estimate = result[2] / pixel_buf.filterWeightSum;
        bool did_isect = std::fabs(result[0]) > 1e-4;

        pixel_buf.xyz[0] = 0;
        pixel_buf.xyz[1] = 0;
        pixel_buf.xyz[2] = 0;

        if ((domain == INTERIOR) == did_isect) {
            FloatToRGB(estimate, result);
            // Hack splatting functionality, since dividing by the number of samples isn't done
            // until the very end for Pixel.xyz
            RGBToXYZ(result, reinterpret_cast<Float*>(pixel_buf.splatXYZ));
            out_file << estimate;
        } else {
            out_file << did_isect;
        }

        if (++i < width) {
            out_file << "\t";
        } else {
            out_file << "\n";
            i = 0;
        }
    }
    out_file << "\n";

    // TODO: write own image reading/writing methods that don't do gamma correction
    camera->film->WriteImage(1);
}

pbrt::Spectrum WoBIntegrator::Li(const pbrt::RayDifferential &r_DO_NOT_USE, const pbrt::Scene &scene,
                                 pbrt::Sampler &sampler, pbrt::MemoryArena &arena, int depth) const
{
    ProfilePhase _(Prof::SamplerIntegratorLi);
    RayDifferential ray(r_DO_NOT_USE);
    Point3f p = ray.o;
    int bounces;
    int start = 0;

    // first, check intersection with the object
    SurfaceInteraction isect_initial;
    auto initial_dir = UniformSampleHemisphere(sampler.Get2D());
    auto did_isect = scene.Intersect(ray, &isect_initial);
    did_isect = did_isect && copysignf(1.0f, isect_initial.p.z) == copysignf(1.0f, p.z);
    if (domain == EXTERIOR && did_isect) {
        did_isect = false;
        ray = isect_initial.SpawnRay(initial_dir);
        p = isect_initial.p;
    } else {
        // intersect with plane z = 0 first, then use that intersection point for the exterior solution
        //  to create effect similar to Sawhney and Crane (2020)
        auto inv_dz = Float(1) / ray.d.z;

        if (isinf(inv_dz)) {
            // camera ray parallel to plane
            return {};
        }

        auto t_ = -p.z * inv_dz;
        p = { p.x + t_ * ray.d.x, p.y + t_ * ray.d.y, 0 };
        ray.o = p;
        ray.d = DomainCoefficient() * initial_dir;

        auto test_ray = Ray(p, r_DO_NOT_USE.d);
        did_isect = did_isect && scene.IntersectP(test_ray);
    }

    if ((domain == INTERIOR) != did_isect) {
        Float res_v[3] = {(Float)did_isect, 0, 0};
        return Spectrum::FromRGB(res_v);
    }

    Float c = DomainCoefficient(); // * Float(scene.OnBoundary(p) ? 2. : 1.);

    Float solution_sample(0.);
    Float pre_solution(0.);
    Float S(1.);

    for (bounces = start; bounces <= maxDepth; ++bounces) {
        // Find next path vertex and accumulate contribution
        VLOG(2) << "Path tracer bounce " << bounces << ", current S = " << S
                << ", solution_sample = " << solution_sample;

        // Intersect _ray_ with scene and store intersection in _isect_
        SurfaceInteraction isect, isect_current;

        size_t m = 0;
        auto o = ray.o;

        // Forward intersection
        while (scene.Intersect(ray, &isect_current)) {
            ray.o = isect_current.p + Float(1e-6) * ray.d;

            m += 1;
            auto u = sampler.Get1D();
            if (Float(m) < Float(1.) / u) {  // equivalent to u < 1/m for reservoir sampling
                isect = isect_current;
            }
        }

        // Backward intersection
        ray.o = o;
        ray.d *= -1;
        ray.tMax = Infinity;

        while (scene.Intersect(ray, &isect_current)) {
            ray.o = isect_current.p + Float(1e-6) * ray.d;

            m += 1;
            auto u = sampler.Get1D();
            if (Float(m) < Float(1.) / u) {
                isect = isect_current;
            }
        }

        // Stop if no intersections found
        if (m == 0) {
            break;
        }

        // Calculate estimate
        // Assume material is MatteMaterial with sigma=0 for LambertianReflection BxDF

        isect.ComputeScatteringFunctions(ray, arena);

        // Sample BSDF to get new path direction
        Vector3f wi = UniformSampleHemisphere(sampler.Get2D());

//        auto ubar = Pi * SpecToFloat(isect.bsdf->f(isect.wo, isect.wo));
        // Manual adjustments for the bunny
        // - 0.25f: the bunny is 8 units wide, since our colourmap is from -1 to 1 we have to scale by 0.25
        // - -0.5f: the bunny's centre is around y = 0.5; this recentres the bunny.
        auto ubar = 0.5f * 0.25f * (isect.p.y - 0.5f);
        // here sign(0) = 1 so that S is unchanged in that case
        S *= Dot(isect.wo, isect.n) < 1e-6 ? Float(m) : -Float(m);
        solution_sample += Float(bounces == maxDepth ? 0.5 : 1) * S * ubar;

        ray = isect.SpawnRay(wi);
    }
    auto res = c * (pre_solution + solution_sample);
    Float res_v[3] = {(Float)did_isect, 128, res};

    return Spectrum::FromRGB(res_v);
}

WoBIntegrator *pbrt_ext::CreateWoBIntegrator(const pbrt::ParamSet &params, std::shared_ptr<pbrt::Sampler> sampler,
                                             std::shared_ptr<const pbrt::Camera> camera)
{
    int maxDepth = params.FindOneInt("maxdepth", 5);
    int np;
    const int *pb = params.FindInt("pixelbounds", &np);
    Bounds2i pixelBounds = camera->film->GetSampleBounds();
    if (pb) {
        if (np != 4)
            Error("Expected four values for \"pixelbounds\" parameter. Got %d.",
                  np);
        else {
            pixelBounds = Intersect(pixelBounds,
                                    Bounds2i{{pb[0], pb[2]}, {pb[1], pb[3]}});
            if (pixelBounds.Area() == 0)
                Error("Degenerate \"pixelbounds\" specified.");
        }
    }
    Float rrThreshold = params.FindOneFloat("rrthreshold", 1.);
    std::string boundaryCond =
        params.FindOneString("boundarycond", "dirichlet");
    std::string domainType =
        params.FindOneString("domaintype", "interior");
    std::string colourmapParam =
        params.FindOneString("colourmap", "turbo");
    return new WoBIntegrator(maxDepth, camera, sampler, pixelBounds,
                             rrThreshold, boundaryCond, domainType, colourmapParam);
}
