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

#define USE_CAMERA 1

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
                             const std::string &colourmapParam)
                             :
#if USE_SAMPLER_INTEGRATOR
    pbrt::SamplerIntegrator(std::move(camera), sampler, pixelBounds),
#else
    camera(std::move(camera)), sampler(std::move(sampler)), pixelBounds(pixelBounds),
#endif
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
    // Render image tiles in parallel

    // Compute number of tiles, _nTiles_, to use for parallel rendering
    Bounds2i sampleBounds = camera->film->GetSampleBounds();
    Vector2i sampleExtent = sampleBounds.Diagonal();
    const int tileSize = 16;
    Point2i nTiles((sampleExtent.x + tileSize - 1) / tileSize,
                   (sampleExtent.y + tileSize - 1) / tileSize);

    // Custom image buffer for the colourmap
    std::unique_ptr<Float[]> img_buf(new Float[3 * sampleBounds.Area()]{0.});

#define PIXEL(p) (p.y * sampleExtent.x + p.x)
#define CHANNEL(p, offset) (3 * PIXEL(p) + offset)

    std::mutex mutex;

#if !USE_SAMPLER_INTEGRATOR
    ProgressReporter reporter(nTiles.x * nTiles.y, "Rendering");
    {
        ParallelFor2D([&](Point2i tile) {
            // Render section of image corresponding to _tile_

            // Allocate _MemoryArena_ for tile
            MemoryArena arena;

            // Get sampler instance for tile
            int seed = tile.y * nTiles.x + tile.x;
            std::unique_ptr<Sampler> tileSampler = sampler->Clone(seed);

            // Compute sample bounds for tile
            int x0 = sampleBounds.pMin.x + tile.x * tileSize;
            int x1 = std::min(x0 + tileSize, sampleBounds.pMax.x);
            int y0 = sampleBounds.pMin.y + tile.y * tileSize;
            int y1 = std::min(y0 + tileSize, sampleBounds.pMax.y);
            Bounds2i tileBounds(Point2i(x0, y0), Point2i(x1, y1));
            LOG(INFO) << "Starting image tile " << tileBounds;

            // Get _FilmTile_ for tile
            std::unique_ptr<FilmTile> filmTile =
                camera->film->GetFilmTile(tileBounds);

            // Loop over pixels in tile to render them
            for (Point2i pixel : tileBounds) {
                {
                    ProfilePhase pp(Prof::StartPixel);
                    tileSampler->StartPixel(pixel);
                }

                // Do this check after the StartPixel() call; this keeps
                // the usage of RNG values from (most) Samplers that use
                // RNGs consistent, which improves reproducability /
                // debugging.
                if (!InsideExclusive(pixel, pixelBounds))
                    continue;

                do {
                    // Initialize _CameraSample_ for current sample
                    CameraSample cameraSample =
                        tileSampler->GetCameraSample(pixel);

                    // Generate camera ray for current sample
                    RayDifferential ray;
#if USE_CAMERA
                    Float rayWeight =
                        camera->GenerateRayDifferential(cameraSample, &ray);
                    ray.ScaleDifferentials(
                        1 / std::sqrt((Float)tileSampler->samplesPerPixel));
#else
                    Float scale = 1;
                    Float px = (Float(pixel.x) + 0.5f) / Float(sampleBounds.Diagonal().x) - 0.5f;
                    Float py = (Float(pixel.y) + 0.5f) / Float(sampleBounds.Diagonal().y) - 0.5f;
                    ray.o = Point3f(scale * px, scale  * py, 0);
                    ray.d = UniformSampleSphere(tileSampler->Get2D());
#endif

                    // Evaluate radiance along camera ray
                    Spectrum L(0.f);
#if USE_CAMERA
                    if (rayWeight > 0) L = Li(ray, scene, *tileSampler, arena);
#else
                    L = Li(ray, scene, *tileSampler, arena);
#endif

                    // Issue warning if unexpected radiance value returned
                    if (L.HasNaNs()) {
                        LOG(ERROR) << StringPrintf(
                            "Not-a-number radiance value returned "
                            "for pixel (%d, %d), sample %d. Setting to black.",
                            pixel.x, pixel.y,
                            (int)tileSampler->CurrentSampleNumber());
                        L = Spectrum(0.);
                    } else if (std::isinf(L.y())) {
                        LOG(ERROR) << StringPrintf(
                            "Infinite luminance value returned "
                            "for pixel (%d, %d), sample %d. Setting to black.",
                            pixel.x, pixel.y,
                            (int)tileSampler->CurrentSampleNumber());
                        L = Spectrum(0.);
                    }
                    VLOG(1) << "Camera sample: " << cameraSample << " -> ray: " <<
                            ray << " -> L = " << L;

                    // Add camera ray's contribution to image
                    // Note: we directly add to the filmTile to bypass the filter pipeline
                    filmTile->GetPixel(pixel).contribSum += L;
                    filmTile->GetPixel(pixel).filterWeightSum += 1.;

                    // Free _MemoryArena_ memory from computing image sample
                    // value
                    arena.Reset();
                } while (tileSampler->StartNextSample());

            }
            LOG(INFO) << "Finished image tile " << tileBounds;

            // Merge image tile into _img_buf_
            {
                ProfilePhase p(Prof::MergeFilmTile);
                LOG(INFO) << "Merging film tile " << tileBounds;
                std::lock_guard<std::mutex> lock(mutex);
                for (Point2i pixel : tileBounds) {
                    // Merge _pixel_ into _Film::pixels_
                    const FilmTilePixel &tilePixel = filmTile->GetPixel(pixel);
                    // still operating in float-only mode, so we only need to add the first channel
                    img_buf[CHANNEL(pixel, 0)] += tilePixel.contribSum[0];
                    img_buf[CHANNEL(pixel, 1)] += tilePixel.filterWeightSum;  // reuse green channel for sample weight
                    img_buf[CHANNEL(pixel, 2)] += tilePixel.contribSum[2];  // use blue channel for object intersection
                }
            }

            reporter.Update();
        }, nTiles);
        reporter.Done();
    }
    LOG(INFO) << "Rendering finished";
#else
    this->pbrt::SamplerIntegrator::Render(scene);

    // Convert results to colourmap
    ParallelFor2D([&](Point2i tile) {
        // Compute sample bounds for tile
        int x0 = sampleBounds.pMin.x + tile.x * tileSize;
        int x1 = std::min(x0 + tileSize, sampleBounds.pMax.x);
        int y0 = sampleBounds.pMin.y + tile.y * tileSize;
        int y1 = std::min(y0 + tileSize, sampleBounds.pMax.y);
        Bounds2i tileBounds(Point2i(x0, y0), Point2i(x1, y1));

        // Get _FilmTile_ for tile
        std::unique_ptr<FilmTile> filmTile =
            camera->film->GetFilmTile(tileBounds);

        // Merge image tile into _img_buf_
        {
            ProfilePhase p(Prof::MergeFilmTile);
            LOG(INFO) << "Merging film tile " << tileBounds;
            std::lock_guard<std::mutex> lock(mutex);
            for (Point2i pixel : tileBounds) {
                // Merge _pixel_ into _Film::pixels_
                const FilmTilePixel &tilePixel = filmTile->GetPixel(pixel);
                // still operating in float-only mode, so we only need to add the first channel
                img_buf[CHANNEL(pixel, 0)] += tilePixel.contribSum[0];
                img_buf[CHANNEL(pixel, 1)] += tilePixel.filterWeightSum;  // reuse green channel for sample weight
                std::cerr << tilePixel.filterWeightSum << std::endl;
            }
        }
    }, nTiles);
#endif

    // Save final image after rendering
    auto centre_pixel = Point2i(sampleBounds.Lerp({0.5f, 0.5f}));
    auto nullspace_correction = cond == NEUMANN ? img_buf[CHANNEL(centre_pixel, 0)] : 0.f;
    std::cerr << "Nullspace correction: " << nullspace_correction / img_buf[CHANNEL(centre_pixel, 1)] << '\n';

    auto width = camera->film->croppedPixelBounds.Diagonal().x;
    auto offset = camera->film->croppedPixelBounds.pMin - sampleBounds.pMin;
#undef PIXEL
#define PIXEL(p) ((p.y - offset.y) * width + p.x - offset.x)

    size_t i = 0;
    std::ofstream out_file(camera->film->filename + ".txt");
    out_file << std::setprecision(3);

    for (auto pixel : camera->film->croppedPixelBounds) {
        Float estimate = img_buf[CHANNEL(pixel, 0)];
        bool did_isect = img_buf[CHANNEL(pixel, 2)] != 0.;
        if (cond == NEUMANN) {
            estimate -= nullspace_correction;
        }
        estimate /= img_buf[CHANNEL(pixel, 1)];
        if ((domain == INTERIOR && did_isect) || (domain == EXTERIOR && !did_isect)) {
            FloatToRGB(estimate, &img_buf[CHANNEL(pixel, 0)]);
            out_file << std::setw(8) << estimate;
        } else {
            img_buf[CHANNEL(pixel, 0)] = 0;
            img_buf[CHANNEL(pixel, 1)] = 0;
            img_buf[CHANNEL(pixel, 2)] = 0;
            out_file << "\t";
        }

        if (++i < sampleExtent.x) {
            out_file << "\t";
        } else {
            out_file << "\n";
            i = 0;
        }
    }
    out_file << "\n";
    // TODO: write own image reading/writing methods that don't do gamma correction
    pbrt::WriteImage(camera->film->filename, img_buf.get(),
                     camera->film->croppedPixelBounds, camera->film->fullResolution);

#undef CHANNEL
#undef PIXEL
}

pbrt::Spectrum WoBIntegrator::Li(const pbrt::RayDifferential &r_DO_NOT_USE, const pbrt::Scene &scene,
                                 pbrt::Sampler &sampler, pbrt::MemoryArena &arena, int depth) const
{
    ProfilePhase _(Prof::SamplerIntegratorLi);
    RayDifferential ray(r_DO_NOT_USE);
    Point3f p = ray.o;
    int bounces;
    int start = 0;

#if USE_CAMERA
    // first, check intersection with the object
    SurfaceInteraction isect_initial;
    auto initial_dir = UniformSampleHemisphere(sampler.Get2D());
    auto did_isect = scene.Intersect(ray, &isect_initial);
    if (domain == EXTERIOR && did_isect && copysignf(1.0f, isect_initial.p.z) == copysignf(1.0f, p.z)) {  // Epsilon
        // check?
        did_isect = true;
        ray = isect_initial.SpawnRay(initial_dir);
        p = isect_initial.p;
    } else {
        did_isect = false;
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

        if (domain == INTERIOR) {
            did_isect = -0.5 <= p.x && p.x <= 0.5 && -0.5 <= p.y && p.y <= 0.5;
        }
    }
#else
    auto did_isect = true;
#endif

    if ((domain == INTERIOR) != did_isect) {
        Float res_v[3] = {0, 0, (Float)did_isect};
        return Spectrum::FromRGB(res_v);
    }

    Float c = DomainCoefficient() * Float(scene.OnBoundary(p) ? 2. : 1.);

    Float solution_sample(0.);
    Float pre_solution(0.);
    Float S(1.);

    for (bounces = start; bounces <= maxDepth; ++bounces) {
        // Find next path vertex and accumulate contribution
        VLOG(2) << "Path tracer bounce " << bounces << ", current S = " << S
                << ", solution_sample = " << solution_sample;

        // Intersect _ray_ with scene and store intersection in _isect_
        SurfaceInteraction isect, isect_current;
        bool foundIntersectionFwd = false, foundIntersectionBckwd = false;

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
            foundIntersectionFwd = true;
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
            foundIntersectionBckwd = true;
        }

        // Stop if no intersections found
        if (!foundIntersectionFwd && !foundIntersectionBckwd) {
            break;
        }

        // Calculate estimate
        // Assume material is MatteMaterial with sigma=0 for LambertianReflection BxDF

        isect.ComputeScatteringFunctions(ray, arena);

        // Sample BSDF to get new path direction
        Vector3f wi = UniformSampleHemisphere(sampler.Get2D());

//        auto ubar = Pi * SpecToFloat(isect.bsdf->f(isect.wo, isect.wo));
        auto ubar = 0.5f * isect.p.y;
        // here sign(0) = 1 so that S is unchanged in that case
        S *= Dot(isect.wo, isect.n) < 1e-6 ? Float(m) : -Float(m);
        solution_sample += Float(bounces == maxDepth ? 0.5 : 1) * S * ubar;

        ray = isect.SpawnRay(wi);
    }
    auto res = c * (pre_solution + solution_sample);
    Float res_v[3] = {res, 16, (Float)did_isect};

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
