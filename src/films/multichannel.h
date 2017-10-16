
/*
    pbrt source code is Copyright(c) 1998-2016
                        Matt Pharr, Greg Humphreys, and Wenzel Jakob.

    This file is part of pbrt.

    Redistribution and use in source and binary forms, with or without
    modification, are permitted provided that the following conditions are
    met:

    - Redistributions of source code must retain the above copyright
      notice, this list of conditions and the following disclaimer.

    - Redistributions in binary form must reproduce the above copyright
      notice, this list of conditions and the following disclaimer in the
      documentation and/or other materials provided with the distribution.

    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
    IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
    TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
    PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
    HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
    SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
    LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
    DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
    THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
    (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
    OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

 */

#if defined(_MSC_VER)
#define NOMINMAX
#pragma once
#endif

#ifndef PBRT_MULTICHANNEL_FILM_H
#define PBRT_MULTICHANNEL_FILM_H

// core/film.h*
#include "film.h"
#include "pbrt.h"
#include "geometry.h"
#include "spectrum.h"
#include "filter.h"
#include "stats.h"
#include "parallel.h"

namespace pbrt {

// MultichannelFilmTilePixel Declarations
    struct MultichannelFilmTilePixel {
        Spectrum contribSum = 0.f;
        Float filterWeightSum = 0.f;
    };

// Film Declarations
    class MultichannelFilm : public Film {
    public:
        // Film Public Methods
        MultichannelFilm(const Point2i &resolution, const Bounds2f &cropWindow,
             std::unique_ptr<Filter> filter, Float diagonal,
             const std::string &filename, Float scale,
             Float maxSampleLuminance = Infinity, const std::string &channelnames, const std::string &renameas);

        void WriteImage(Float splatScale = 1);

//    private:
//        // Film Private Data
//        struct Pixel {
//            Pixel() { xyz[0] = xyz[1] = xyz[2] = filterWeightSum = 0; }
//            Float xyz[3];
//            Float filterWeightSum;
//            AtomicFloat splatXYZ[3];
//            Float pad;
//        };
//        std::unique_ptr<Pixel[]> pixels;
//        static PBRT_CONSTEXPR int filterTableWidth = 16;
//        Float filterTable[filterTableWidth * filterTableWidth];
//        std::mutex mutex;
//        const Float scale;
//        const Float maxSampleLuminance;
//
//        // Film Private Methods
//        Pixel &GetPixel(const Point2i &p) {
//            CHECK(InsideExclusive(p, croppedPixelBounds));
//            int width = croppedPixelBounds.pMax.x - croppedPixelBounds.pMin.x;
//            int offset = (p.x - croppedPixelBounds.pMin.x) +
//                         (p.y - croppedPixelBounds.pMin.y) * width;
//            return pixels[offset];
//        }
    };

    class MultichannelFilmTile : public FilmTile{
    public:
        // FilmTile Public Methods
        MultichannelFilmTile(const Bounds2i &pixelBounds, const Vector2f &filterRadius,
                 const Float *filterTable, int filterTableSize,
                 Float maxSampleLuminance)
                : pixelBounds(pixelBounds),
                  filterRadius(filterRadius),
                  invFilterRadius(1 / filterRadius.x, 1 / filterRadius.y),
                  filterTable(filterTable),
                  filterTableSize(filterTableSize),
                  maxSampleLuminance(maxSampleLuminance) {
            pixels = std::vector<MultichannelFilmTilePixel>(std::max(0, pixelBounds.Area()));
        }

    private:
        // FilmTile Private Data
        const Bounds2i pixelBounds;
        const Vector2f filterRadius, invFilterRadius;
        const Float *filterTable;
        const int filterTableSize;
        std::vector<MultichannelFilmTilePixel> pixels;
        const Float maxSampleLuminance;
        friend class MultichannelFilm;
    };

    MultichannelFilm *CreateFilm(const ParamSet &params, std::unique_ptr<Filter> filter);

}  // namespace pbrt

#endif  // PBRT_CORE_FILM_H
