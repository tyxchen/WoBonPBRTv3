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

WoBIntegrator::WoBIntegrator(int maxDepth, std::shared_ptr<const Camera> camera,
                             std::shared_ptr<Sampler> sampler,
                             const pbrt::Bounds2i &pixelBounds, pbrt::Float rrThreshold,
                             const std::string &boundaryCond,
                             const std::string &domainType,
                             const std::string &colourmapParam) :
    pbrt::SamplerIntegrator(std::move(camera), sampler, pixelBounds),
    maxDepth(maxDepth), rrThreshold(rrThreshold)
{
    if (boundaryCond == "dirichlet") {
        cond = DIRICHLET;
    } else {
        throw std::logic_error("Bad boundary condition given.");
    }

    if (domainType == "interior") {
        domain = INTERIOR;
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
    std::ofstream mask_file(camera->film->filename + ".mask");
    out_file << std::setprecision(3);

    for (auto pixel : camera->film->croppedPixelBounds) {
        auto &pixel_buf = camera->film->GetPixel(pixel);
        Float result[3];
        XYZToRGB(pixel_buf.xyz, result);
        auto estimate = result[1] / pixel_buf.filterWeightSum;
        bool in_interior = std::fabs(result[0]) != 0.;

        pixel_buf.xyz[0] = 0;
        pixel_buf.xyz[1] = 0;
        pixel_buf.xyz[2] = 0;

        if ((domain == INTERIOR) == in_interior) {
            FloatToRGB(estimate, result);
            // Hack splatting functionality, since dividing by the number of samples isn't done
            // until the very end for Pixel.xyz
            RGBToXYZ(result, reinterpret_cast<Float*>(pixel_buf.splatXYZ));
            out_file << estimate;
        } else {
            pixel_buf.splatXYZ[0] = 1;
            pixel_buf.splatXYZ[1] = 1;
            pixel_buf.splatXYZ[2] = 1;
            out_file << 0;
        }

        mask_file << in_interior;

        if (++i < width) {
            out_file << ",";
            mask_file << ",";
        } else {
            out_file << "\n";
            mask_file << "\n";
            i = 0;
        }
    }
    out_file << "\n";
    mask_file << "\n";

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
    bool in_interior = CheckInteriorPoint(scene, sampler, p, ray);

    if ((domain == INTERIOR) != in_interior) {
        Float res_v[3] = { (Float)in_interior, 0, 0};
        return Spectrum::FromRGB(res_v);
    }

    Float solution_sample(0.);
    Float pre_solution(0.);
    Float S(1.);
    SurfaceInteraction isect;

    for (bounces = start; bounces <= maxDepth; ++bounces) {
        // Find next path vertex and accumulate contribution
        VLOG(2) << "Path tracer bounce " << bounces << ", current S = " << S
                << ", solution_sample = " << solution_sample;

        // Intersect _ray_ with scene and store intersection in _isect_
        RayDifferential sampled_ray;

        size_t m = 0;
        auto prev_isect = isect;
        auto d = ray.d;

        // Forward intersection
        if (bounces > start) {
            ray = prev_isect.SpawnRay(d);
        }
        Intersect(scene, ray, sampler, isect, sampled_ray, m);

        // Backward intersection
        if (bounces > start) {
            ray = prev_isect.SpawnRay(-d);
        } else {
            ray = Ray(p, -d);
        }
        Intersect(scene, ray, sampler, isect, sampled_ray, m);

        // Stop if no intersections found
        if (m == 0) {
            break;
        }

        // Calculate estimate
        isect.ComputeScatteringFunctions(sampled_ray, arena);

        // Sample BSDF to get new path direction
        Vector3f wi = UniformSampleHemisphere(sampler.Get2D());

        auto ubar = isect.p.y;
        // here sign(0) = 1 so that S is unchanged in that case
        S *= Dot(sampled_ray.d, isect.n) < 1e-6 ? Float(m) : -Float(m);
        solution_sample += Float(bounces == maxDepth ? 0.5 : 1) * S * ubar;

        ray.d = wi;
        p = isect.p;
    }

    auto res = DomainCoefficient() * (pre_solution + solution_sample);
    // [ mask, numerical solution, dummy value ]
    Float res_v[3] = { (Float)in_interior, res, 0};

    return Spectrum::FromRGB(res_v);
}

// Any-hit intersection routine
void WoBIntegrator::Intersect(const Scene &scene, RayDifferential ray, Sampler &sampler, SurfaceInteraction &isect,
                              RayDifferential &sampled_ray, size_t &m) const
{
    SurfaceInteraction isect_current;
    while (scene.Intersect(ray, &isect_current)) {
        ray = isect_current.SpawnRay(ray.d);

        m += 1;
        auto u = sampler.Get1D();
        if (u < Float(1.) / Float(m)) {  // equivalent to u < 1/m for reservoir sampling
            isect = isect_current;
            sampled_ray = ray;
        }
    }
}

// Check if p is in the interior of the domain and generate the initial sample ray
bool WoBIntegrator::CheckInteriorPoint(const Scene &scene, Sampler &sampler, Point3f &p, RayDifferential &ray) const
{
    SurfaceInteraction isect_initial;
    bool in_interior = scene.Intersect(ray, &isect_initial);
    auto initial_dir = UniformSampleHemisphere(sampler.Get2D());
    in_interior = in_interior && copysignf(1.0f, isect_initial.p.z) == copysignf(1.0f, p.z);
    // intersect with plane z = 0 first, then use that intersection point for the exterior solution
    //  to create effect similar to Sawhney and Crane (2020)
    auto inv_dz = Float(1) / ray.d.z;

    auto t_ = -p.z * inv_dz;
    p = Point3f(p.x + t_ * ray.d.x, p.y + t_ * ray.d.y, 0);
    ray = Ray(p, DomainCoefficient() * initial_dir);

    if (in_interior) {
        auto test_ray = Ray(p, Normalize(Vector3f(p)));
        size_t num_hits = 0;
        SurfaceInteraction interior_isect;

        while (scene.Intersect(test_ray, &interior_isect)) {
            test_ray = interior_isect.SpawnRay(test_ray.d);
            ++num_hits;
        }

        in_interior = num_hits % 2 == 1;
    }

    return in_interior;
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
    std::cout << "Recursion Depth: " << maxDepth << std::endl;
    return new WoBIntegrator(maxDepth, camera, sampler, pixelBounds,
                             rrThreshold, boundaryCond, domainType, colourmapParam);
}
