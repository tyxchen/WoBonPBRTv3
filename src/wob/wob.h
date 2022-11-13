//
// Created by Terry Chen on 2022-10-21.
//

#ifndef PBRT_V3_WOB_H
#define PBRT_V3_WOB_H

#include "pbrt.h"
#include "integrator.h"
#include "lightdistrib.h"

namespace pbrt_ext
{

class WoBIntegrator : public pbrt::Integrator
{
    enum BoundaryCond {
        DIRICHLET,
        NEUMANN,
        ROBIN,
    };
public:
    WoBIntegrator(int maxDepth,
                  std::shared_ptr<const pbrt::Camera> camera,
                  std::shared_ptr<pbrt::Sampler> sampler,
                  const pbrt::Bounds2i &pixelBounds, pbrt::Float rrThreshold = 1,
                  const std::string &boundaryCond = "dirichlet");

    pbrt::Float Li(const pbrt::RayDifferential &ray, const pbrt::Scene &scene, pbrt::Sampler &sampler,
                   pbrt::MemoryArena &arena, int depth = 0) const;

    void Render(const pbrt::Scene &scene) override;

private:
    std::shared_ptr<pbrt::Sampler> sampler;
    const pbrt::Bounds2i pixelBounds;
    std::shared_ptr<const pbrt::Camera> camera;
    int maxDepth;
    BoundaryCond cond;
};

WoBIntegrator *CreateWoBIntegrator(const pbrt::ParamSet &params,
                                   std::shared_ptr<pbrt::Sampler> sampler,
                                   std::shared_ptr<const pbrt::Camera> camera);

} // pbrt_ext

#endif //PBRT_V3_WOB_H
