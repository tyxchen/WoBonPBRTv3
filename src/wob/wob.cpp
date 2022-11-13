//
// Created by Terry Chen on 2022-10-21.
//

#include "core/camera.h"
#include "core/film.h"
#include "core/progressreporter.h"
#include "core/scene.h"
#include "core/paramset.h"
#include "wob/wob.h"

using namespace pbrt;
using namespace pbrt_ext;

inline Float G(const Point3f &x, const Point3f &y) {
    return Inv4Pi / (x - y).Length();
}

inline Float H(const Point3f &x, const Point3f &y, const Normal3f &n_y) {
    auto r = (y - x).Length();
    return -Inv4Pi * Dot(y - x, n_y) / (r * r * r);
}

inline Float SpecToFloat(const Spectrum &s) {
    return Pi * (s[0] + -s[2]);
}

inline Spectrum FloatToSpec(Float f) {
    constexpr
    #include <wob/turbo_colormap.c>

    auto clamped_f = std::min(std::max(f, Float(-1.)), Float(1.));
    clamped_f = Float(0.5) * (clamped_f + Float(1.));
    return Spectrum::FromRGB(turbo_srgb_floats[(int)std::floor(clamped_f * 255)]);
}

WoBIntegrator::WoBIntegrator(int maxDepth, std::shared_ptr<const Camera> camera,
                             std::shared_ptr<Sampler> sampler,
                             const pbrt::Bounds2i &pixelBounds, pbrt::Float rrThreshold,
                             const std::string &boundaryCond)
: sampler(sampler), pixelBounds(pixelBounds), camera(camera), maxDepth(5)
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
                    Float px = Float(pixel.x) / Float(sampleBounds.Diagonal().x);
                    Float py = Float(pixel.y) / Float(sampleBounds.Diagonal().y);
                    ray.o = Point3f(px, py, 0.);
                    ray.d = UniformSampleHemisphere(tileSampler->Get2D());

//                    ray.ScaleDifferentials(
//                        1 / std::sqrt((Float)tileSampler->samplesPerPixel));

                    // Evaluate radiance along camera ray
                    Float L = Li(ray, scene, *tileSampler, arena);

                    // Issue warning if unexpected radiance value returned
                    if (std::isnan(L)) {
                        LOG(ERROR) << StringPrintf(
                            "Not-a-number radiance value returned "
                            "for pixel (%d, %d), sample %d. Setting to black.",
                            pixel.x, pixel.y,
                            (int)tileSampler->CurrentSampleNumber());
                        L = Float(0.);
                    } else if (std::isinf(L)) {
                        LOG(ERROR) << StringPrintf(
                            "Infinite luminance value returned "
                            "for pixel (%d, %d), sample %d. Setting to black.",
                            pixel.x, pixel.y,
                            (int)tileSampler->CurrentSampleNumber());
                        L = Float(0.);
                    }
                    VLOG(1) << "Camera sample: " << cameraSample << " -> ray: " <<
                            ray << " -> L = " << L;

                    // Add camera ray's contribution to image
                    filmTile->AddSample(cameraSample.pFilm, Spectrum(L), 1.);

                    // Free _MemoryArena_ memory from computing image sample
                    // value
                    arena.Reset();
                } while (tileSampler->StartNextSample());

                filmTile->GetPixel(pixel).contribSum = FloatToSpec(filmTile->GetPixel(pixel).contribSum[0] /
                    Float(tileSampler->samplesPerPixel));
            }
            LOG(INFO) << "Finished image tile " << tileBounds;

            // Merge image tile into _Film_
            camera->film->MergeFilmTile(std::move(filmTile));
            reporter.Update();
        }, nTiles);
        reporter.Done();
    }
    LOG(INFO) << "Rendering finished";

    // Save final image after rendering
    camera->film->WriteImage();
}

Float WoBIntegrator::Li(const pbrt::RayDifferential &r, const pbrt::Scene &scene, pbrt::Sampler &sampler,
                        pbrt::MemoryArena &arena, int depth) const
{
    ProfilePhase _(Prof::SamplerIntegratorLi);
    Point3f p = r.o;
    RayDifferential ray(r);
    int bounces;

    Float c = scene.OnBoundary(p) ? 2. : 1.;

    Float solution_sample(0.);
    Float pre_solution(0.);
    Float S(1.);

    if (cond == NEUMANN) {
        // calculate pre-solution for NEE
        Float pdf_sample;
        auto vhit = scene.Sample(sampler.Get2D(), &pdf_sample);
        vhit.ComputeScatteringFunctions(r, arena);
        auto qbar = SpecToFloat(vhit.bsdf->f(vhit.wo, vhit.wo));
        pre_solution += G(p, vhit.p) * qbar / pdf_sample;
    }

    for (bounces = 0; bounces <= maxDepth; ++bounces) {
        // Find next path vertex and accumulate contribution
        VLOG(2) << "Path tracer bounce " << bounces << ", current S = " << S
                << ", solution_sample = " << solution_sample;

        // Intersect _ray_ with scene and store intersection in _isect_
        SurfaceInteraction isect, isect_current;
        bool foundIntersectionFwd = false, foundIntersectionBckwd = false;

        size_t m = 0;
        Float pdf_isect;
        auto o = ray.o;

        while (scene.Intersect(ray, &isect_current)) {
            ray.o = isect_current.p + Float(1e-6) * ray.d;

            pdf_isect = Float(1) / Float(++m);
            auto u = sampler.Get1D();
            if (u < pdf_isect) {
                isect = isect_current;
            }
            foundIntersectionFwd = true;
        }

        ray.o = o;
        ray.d *= -1;
        ray.tMax = Infinity;

        while (scene.Intersect(ray, &isect_current)) {
            ray.o = isect_current.p + Float(1e-6) * ray.d;

            pdf_isect = Float(1) / Float(++m);
            auto u = sampler.Get1D();
            if (u < pdf_isect) {
                isect = isect_current;
            }
            foundIntersectionBckwd = true;
        }

        if (!foundIntersectionFwd && !foundIntersectionBckwd) {
            break;
        }

        // Calculate estimate
        // Assume material is MatteMaterial with sigma=0 for LambertianReflection BxDF

        switch (cond) {
        case DIRICHLET: {
            isect.ComputeScatteringFunctions(r, arena);
            auto ubar = SpecToFloat(isect.bsdf->f(isect.wo, isect.wo));
            // here sign(0) = -1 so that S is unchanged in that case
            S *= -Float(Dot(isect.wo, isect.n) < -1e-6 ? m : -m);
            solution_sample += Float(bounces == maxDepth ? 0.5 : 1) * S * ubar;
            break;
        }
        case NEUMANN: {
            Float pdf_sample;
            auto vhit = scene.Sample(sampler.Get2D(), &pdf_sample);
            vhit.ComputeScatteringFunctions(r, arena);
            auto qbar = SpecToFloat(vhit.bsdf->f(vhit.wo, vhit.wo) * Pi);
            S *= Float(Dot(isect.wo, isect.n) < 1e-6 ? m : -m);
            solution_sample += S * G(isect.p, vhit.p) * qbar / pdf_sample;
            break;
        }
        case ROBIN: {
            // TODO: implement
            break;
        }
        }

        // Sample BSDF to get new path direction
        auto wi = UniformSampleHemisphere(sampler.Get2D());
        ray = isect.SpawnRay(wi);

    }
    return c * (pre_solution + solution_sample);
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
    return new WoBIntegrator(maxDepth, camera, sampler, pixelBounds,
                             rrThreshold, boundaryCond);
}
