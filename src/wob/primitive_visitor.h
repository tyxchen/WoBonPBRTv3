//
// Created by Terry Chen on 2022-11-11.
//

#ifndef PBRT_V3_PRIMITIVE_VISITOR_H
#define PBRT_V3_PRIMITIVE_VISITOR_H

#include "core/pbrt.h"

namespace pbrt_ext
{
template <typename T>
struct PrimitiveVisitor {
    inline void build(const std::vector<std::shared_ptr<pbrt::Primitive>> *prims) {
        primitives = prims;
        // Calculate cumulative area
        cumulativeArea.reserve(primitives->size() + 1);
        cumulativeArea.push_back(0);
        for (const auto &prim : *primitives) {
            cumulativeArea.push_back(cumulativeArea.back() + prim->Area());
        }
    }

    inline pbrt::SurfaceInteraction sample(const pbrt::Point2f &u, pbrt::Float *pdf) const {
        size_t n = primitives->size();
        auto v = cumulativeArea[n] * u.x;
        size_t lower = 0, upper = n;
        while (lower < upper - 1) {
            auto mid = (upper + lower) / 2;
            if (v < cumulativeArea[mid]) {
                upper = mid;
            } else {
                lower = mid;
            }
        }
        const auto f = (*primitives)[lower];
        v = (v - cumulativeArea[lower]) / (cumulativeArea[upper] - cumulativeArea[lower]);
        *pdf = pbrt::Float(1.) / cumulativeArea[n];
        pbrt::Float dummy;
        return f->Sample({v, u.y}, &dummy);
    }

    inline pbrt::Float winding_number(const pbrt::Point3f &p) const {
        // Winding number calculation from Jacobson, Kavan, and Sorkine-Hornung 2013
        pbrt::Float winding = 0;
        for (const auto &prim : *primitives) {
            winding += prim->SolidAngle(p);
        }
        return winding * pbrt::Inv4Pi;
    }

    std::vector<pbrt::Float> cumulativeArea;
    const std::vector<std::shared_ptr<pbrt::Primitive>> *primitives = nullptr;
};
}

#endif //PBRT_V3_PRIMITIVE_VISITOR_H
