//
// Created by terry on 21/01/23.
//

#ifndef PBRT_V3_STRIPETEXTURE_H
#define PBRT_V3_STRIPETEXTURE_H

#include "pbrt.h"

class StripeTexture : public pbrt::Texture<pbrt::Spectrum>
{
public:
    StripeTexture(const pbrt::Spectrum &K0, const pbrt::Spectrum &K1) : K0(K0), K1(K1) {}
    pbrt::Spectrum Evaluate(const pbrt::SurfaceInteraction &si) const override {
        auto yint = (int)floor(si.p.y);
        auto is_stripe = (yint + 8) % 2 == 1;
        return is_stripe ? K1 : K0;
    }
private:
    pbrt::Spectrum K0, K1;
};

StripeTexture *CreateStripeSpectrumTexture(const pbrt::Transform &tex2world, const pbrt::TextureParams &tp) {
    pbrt::Spectrum K0 = tp.FindSpectrum("K0", pbrt::Spectrum());
    pbrt::Spectrum K1 = tp.FindSpectrum("K1", pbrt::Spectrum(1.f));
    return new StripeTexture(K0, K1);
}

#endif //PBRT_V3_STRIPETEXTURE_H
