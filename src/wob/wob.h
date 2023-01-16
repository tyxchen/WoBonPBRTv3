//
// Created by Terry Chen on 2022-10-21.
//

#ifndef PBRT_V3_WOB_H
#define PBRT_V3_WOB_H

#include "pbrt.h"
#include "integrator.h"
#include "lightdistrib.h"

#include <tinycolormap.hpp>

#define USE_SAMPLER_INTEGRATOR 0

namespace pbrt_ext
{

#if USE_SAMPLER_INTEGRATOR
class WoBIntegrator : public pbrt::SamplerIntegrator
#else
class WoBIntegrator : public pbrt::Integrator
#endif
{
    enum BoundaryCond {
        DIRICHLET,
        NEUMANN,
        ROBIN,
    };
    enum DomainType {
        INTERIOR,
        EXTERIOR,
    };
public:
    WoBIntegrator(int maxDepth,
                  std::shared_ptr<const pbrt::Camera> camera,
                  std::shared_ptr<pbrt::Sampler> sampler,
                  const pbrt::Bounds2i &pixelBounds, pbrt::Float rrThreshold = 1,
                  const std::string &boundaryCond = "dirichlet",
                  const std::string &domainType = "interior",
                  const std::string &colourmapParam = "turbo");

    pbrt::Spectrum Li(const pbrt::RayDifferential &ray, const pbrt::Scene &scene, pbrt::Sampler &sampler,
                      pbrt::MemoryArena &arena, int depth = 0) const;

    void Render(const pbrt::Scene &scene) override;
private:
#if !USE_SAMPLER_INTEGRATOR
    std::shared_ptr<pbrt::Sampler> sampler;
    const pbrt::Bounds2i pixelBounds;
    std::shared_ptr<const pbrt::Camera> camera;
#endif
    const int maxDepth;
    const pbrt::Float rrThreshold;
    BoundaryCond cond;
    DomainType domain;
    tinycolormap::ColormapType colourmap;

    inline pbrt::Float SpecToFloat(const pbrt::Spectrum &s) const {
        return s[0] + -s[2];
    }

    inline void FloatToRGB(pbrt::Float f, pbrt::Float *out) const {
        f = pbrt::Float(0.5) * (f + pbrt::Float(1.));
        auto mapped_colour = tinycolormap::GetColor(f, colourmap);
        out[0] = static_cast<pbrt::Float>(mapped_colour.data[0]);
        out[1] = static_cast<pbrt::Float>(mapped_colour.data[1]);
        out[2] = static_cast<pbrt::Float>(mapped_colour.data[2]);
    }

    inline pbrt::Float DomainCoefficient() const {
        return domain == INTERIOR ? 1 : -1;
    }
};

WoBIntegrator *CreateWoBIntegrator(const pbrt::ParamSet &params,
                                   std::shared_ptr<pbrt::Sampler> sampler,
                                   std::shared_ptr<const pbrt::Camera> camera);

} // pbrt_ext

#endif //PBRT_V3_WOB_H
