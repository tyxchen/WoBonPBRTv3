//
// Created by Terry Chen on 2022-11-11.
//

#ifndef PBRT_V3_PRIMITIVE_VISITOR_H
#define PBRT_V3_PRIMITIVE_VISITOR_H

#include "core/pbrt.h"

namespace pbrt_ext
{
template <typename T>
class PrimitiveVisitor {
    friend T;

    inline static void build(T * const host) {
        // Calculate cumulative area
        host->cumulativeArea.reserve(host->primitives.size() + 1);
        host->cumulativeArea.push_back(0);
        for (const auto &prim : host->primitives) {
            host->cumulativeArea.push_back(host->cumulativeArea.back() + prim->Area());
        }
    }

    inline static pbrt::SurfaceInteraction sample(const T * const host, const pbrt::Point2f &u, pbrt::Float *pdf) {
        size_t n = host->primitives.size();
        auto v = host->cumulativeArea[n] * u.x;
        size_t lower = 0, upper = n;
        while (lower < upper - 1) {
            auto mid = (upper + lower) / 2;
            if (v < host->cumulativeArea[mid]) {
                upper = mid;
            } else {
                lower = mid;
            }
        }
        const auto f = host->primitives[lower];
        v = (v - host->cumulativeArea[lower]) / (host->cumulativeArea[upper] - host->cumulativeArea[lower]);
        *pdf = pbrt::Float(1.) / host->cumulativeArea[n];
        pbrt::Float dummy;
        return f->Sample({v, u.y}, &dummy);
    }

    inline static pbrt::Float winding_number(const T * const host, const pbrt::Point3f &p) {
        // Winding number calculation from Jacobson, Kavan, and Sorkine-Hornung 2013
        pbrt::Float winding = 0;
        for (const auto &prim : host->primitives) {
            winding += prim->SolidAngle(p);
        }
        return winding * pbrt::Inv4Pi;
    }
};
}

#endif //PBRT_V3_PRIMITIVE_VISITOR_H
