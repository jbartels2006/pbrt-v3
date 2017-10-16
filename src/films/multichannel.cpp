
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


// core/film.cpp*
#include "multichannel.h"
#include "paramset.h"
#include "imageio.h"
#include "stats.h"

namespace pbrt {

    STAT_MEMORY_COUNTER("Memory/Film pixels", filmPixelMemory);

// Film Method Definitions
    MultichannelFilm::MultichannelFilm(const Point2i &resolution, const Bounds2f &cropWindow,
               std::unique_ptr<Filter> filt, Float diagonal,
               const std::string &filename, Float scale, Float maxSampleLuminance, const std::string &channelnames, const std::string &renameas)
            : fullResolution(resolution),
              diagonal(diagonal * .001),
              filter(std::move(filt)),
              filename(filename),
              scale(scale),
              maxSampleLuminance(maxSampleLuminance) {
        // Compute film image bounds
        croppedPixelBounds =
                Bounds2i(Point2i(std::ceil(fullResolution.x * cropWindow.pMin.x),
                                 std::ceil(fullResolution.y * cropWindow.pMin.y)),
                         Point2i(std::ceil(fullResolution.x * cropWindow.pMax.x),
                                 std::ceil(fullResolution.y * cropWindow.pMax.y)));
        LOG(INFO) << "Created film with full resolution " << resolution <<
                  ". Crop window of " << cropWindow << " -> croppedPixelBounds " <<
                  croppedPixelBounds;

        // Allocate film image storage
        pixels = std::unique_ptr<Pixel[]>(new Pixel[croppedPixelBounds.Area()]);
        filmPixelMemory += croppedPixelBounds.Area() * sizeof(Pixel);

        // Precompute filter weight table
        int offset = 0;
        for (int y = 0; y < filterTableWidth; ++y) {
            for (int x = 0; x < filterTableWidth; ++x, ++offset) {
                Point2f p;
                p.x = (x + 0.5f) * filter->radius.x / filterTableWidth;
                p.y = (y + 0.5f) * filter->radius.y / filterTableWidth;
                filterTable[offset] = filter->Evaluate(p);
            }
        }
    }


    // Look up the channel name and return the appropriate value.  If the channel name
// is not found, we show a warning and quit.
    float MultichannelFilm::ComputeVal(const std::string &channelname, const Ray &ray) const
    {
        std::string s = channelname;

        transform (s.begin(), s.end(), s.begin(), tolower);
        if (s == "z") return (float) ray.tMax;
        if (s == "wx") return ray.p.x;
        if (s == "wy") return ray.p.y;
        if (s == "wz") return ray.p.z;
        if (s == "nx") return ray.n.x;
        if (s == "ny") return ray.n.y;
        if (s == "nz") return ray.n.z;
        if (s == "ou") return ray.uv.x;
        if (s == "ov") return ray.uv.y;

        // The specified channel does not exist.  Exit and show a warning.
        std::cout << "\n";
        Warning("Channel '%s' does not exist.  Exiting...", channelname.c_str());

        exit(1);
    }

    void MultichannelFilm::UpdateChannel(const string &channelname, int index, const RayDifferential &ray, int x, int y, float weight)
    {
        float val = ComputeVal(channelname, ray);

        // We will arbitrarily choose the last sample that we looked at.  (Ideally we would
        // keep a histogram of all values and choose the mode...)
        (*nonRgbaChannelData)((x * numNonRgbaChannels) + (index), y) = val;
    }


    Imf::PixelType MultichannelFilm::GetVarType(const std::string &channelname)
    {
        if (channelname == "Z") return Imf::FLOAT;

        return Imf::HALF;
    }

    void MultichannelFilm::WriteImage(Float splatScale) {
        // Convert image to RGB and compute final pixel values
        LOG(INFO) <<
                  "Converting image to RGB and computing final weighted pixel values";
        std::unique_ptr<Float[]> rgb(new Float[3 * croppedPixelBounds.Area()]);
        int offset = 0;
        for (Point2i p : croppedPixelBounds) {
            // Convert pixel XYZ color to RGB
            Pixel &pixel = GetPixel(p);
            XYZToRGB(pixel.xyz, &rgb[3 * offset]);

            // Normalize pixel with weight sum
            Float filterWeightSum = pixel.filterWeightSum;
            if (filterWeightSum != 0) {
                Float invWt = (Float)1 / filterWeightSum;
                rgb[3 * offset] = std::max((Float)0, rgb[3 * offset] * invWt);
                rgb[3 * offset + 1] =
                        std::max((Float)0, rgb[3 * offset + 1] * invWt);
                rgb[3 * offset + 2] =
                        std::max((Float)0, rgb[3 * offset + 2] * invWt);
            }

            // Add splat value at pixel
            Float splatRGB[3];
            Float splatXYZ[3] = {pixel.splatXYZ[0], pixel.splatXYZ[1],
                                 pixel.splatXYZ[2]};
            XYZToRGB(splatXYZ, splatRGB);
            rgb[3 * offset] += splatScale * splatRGB[0];
            rgb[3 * offset + 1] += splatScale * splatRGB[1];
            rgb[3 * offset + 2] += splatScale * splatRGB[2];

            // Scale pixel value by _scale_
            rgb[3 * offset] *= scale;
            rgb[3 * offset + 1] *= scale;
            rgb[3 * offset + 2] *= scale;
            ++offset;

            //TODO: Add code to get multichannel stuff here?
            // Also write non-RGBA channels
            RayDifferential &ray = GetPixelRay(p);
            for (int i = 0; i < numNonRgbaChannels; i++)
            {
                UpdateChannel(nonRgbaChannelNames[i], nonRgbaChannelIndex[i], ray, p.x, p.y, filterWt);
            }
        }

        // Write RGB image
        LOG(INFO) << "Writing image " << filename << " with bounds " <<
                  croppedPixelBounds;
        pbrt::WriteImage(filename, &rgb[0], croppedPixelBounds, fullResolution); //TODO: change to writeextendedchannel image
    }

    void MultichannelFilm::WriteExtendedChannelImage(const std::string &name, const std::vector<std::string> &channelnames,
                                                        const std::vector<Imf::PixelType> &types, float *pixels,
                                                        int xRes, int yRes,
                                                        int totalXRes, int totalYRes,
                                                        int xOffset, int yOffset)
    {
        Imf::Header header(totalXRes, totalYRes);
        Box2i dataWindow(V2i(xOffset, yOffset), V2i(xOffset + xRes - 1, yOffset + yRes - 1));
        header.dataWindow() = dataWindow;

        int numchannels = channelnames.size();
        for (int i = 0; i < numchannels; i++)
        {
            header.channels().insert(channelnames[i].c_str(), Imf::Channel (types[i]));
        }

        // Copy all data for each channel into a single array and then
        // add that into the framebuffer
        Imf::FrameBuffer fb;

        half **hptrs = new half*[numchannels];
        unsigned int **iptrs = new unsigned int *[numchannels];
        float **fptrs = new float*[numchannels];

        for (int i = 0; i < numchannels; i++)
        {
            hptrs[i] = NULL;
            fptrs[i] = NULL;
            iptrs[i] = NULL;
        }

        for (int i = 0; i < numchannels; i++)
        {
            if (types[i] == Imf::HALF)
            {
                half *h = new half[xRes * yRes];
                hptrs[i] = h;

                for (int j = 0; j < xRes * yRes; j++)
                {
                    h[j] = (half) pixels[j * numchannels + i];
                }

                h -= (xOffset + yOffset * xRes);

                fb.insert(channelnames[i].c_str(), Imf::Slice(types[i], (char *)h, sizeof(half),
                                                              xRes*sizeof(half)));

            } else if ((types[i]) == Imf::FLOAT)
            {
                float *f = new float[xRes * yRes];
                fptrs[i] = f;

                for (int j = 0; j < xRes * yRes; j++)
                {
                    f[j] = pixels[j * numchannels + i];
                }

                f -= (xOffset + yOffset * xRes);

                fb.insert(channelnames[i].c_str(), Imf::Slice(types[i], (char *)f, sizeof(float),
                                                              xRes*sizeof(float)));

            } else
            {
                unsigned int *ui = new unsigned int[xRes * yRes];
                iptrs[i] = ui;

                for (int j = 0; j < xRes * yRes; j++)
                {
                    ui[j] = (unsigned int) pixels[j * numchannels + i];
                }

                ui -= (xOffset + yOffset * xRes);

                fb.insert(channelnames[i].c_str(), Imf::Slice(types[i], (char *)ui, sizeof(unsigned int),
                                                              xRes*sizeof(unsigned int)));

            }
        }

        Imf::OutputFile file(name.c_str(), header);

        file.setFrameBuffer(fb);
        try {
            file.writePixels(yRes);
        }
        catch (const std::exception &e) {
            Error("Unable to write image file \"%s\": %s", name.c_str(),
                  e.what());
        }
        // Free the memory that was allocated for the slices
        for (int i = 0; i < numchannels; i++)
        {
            delete[] fptrs[i];
            delete[] iptrs[i];
            delete[] hptrs[i];
        }
        delete[] fptrs;
        delete[] iptrs;
        delete[] hptrs;
    }

    MultichannelFilm *CreateFilm(const ParamSet &params, std::unique_ptr<Filter> filter) {
        // Intentionally use FindOneString() rather than FindOneFilename() here
        // so that the rendered image is left in the working directory, rather
        // than the directory the scene file lives in.
        std::string filename = params.FindOneString("filename", "");
        if (PbrtOptions.imageFile != "") {
            if (filename != "") {
                Warning(
                        "Output filename supplied on command line, \"%s\", ignored "
                                "due to filename provided in scene description file, \"%s\".",
                        PbrtOptions.imageFile.c_str(), filename.c_str());
            } else
                filename = PbrtOptions.imageFile;
        }
        if (filename == "") filename = "pbrt.exr";

        int xres = params.FindOneInt("xresolution", 1280);
        int yres = params.FindOneInt("yresolution", 720);
        if (PbrtOptions.quickRender) xres = std::max(1, xres / 4);
        if (PbrtOptions.quickRender) yres = std::max(1, yres / 4);
        Bounds2f crop(Point2f(0, 0), Point2f(1, 1));
        int cwi;
        const Float *cr = params.FindFloat("cropwindow", &cwi);
        if (cr && cwi == 4) {
            crop.pMin.x = Clamp(std::min(cr[0], cr[1]), 0.f, 1.f);
            crop.pMax.x = Clamp(std::max(cr[0], cr[1]), 0.f, 1.f);
            crop.pMin.y = Clamp(std::min(cr[2], cr[3]), 0.f, 1.f);
            crop.pMax.y = Clamp(std::max(cr[2], cr[3]), 0.f, 1.f);
        } else if (cr)
            Error("%d values supplied for \"cropwindow\". Expected 4.", cwi);

        Float scale = params.FindOneFloat("scale", 1.);
        Float diagonal = params.FindOneFloat("diagonal", 35.);
        Float maxSampleLuminance = params.FindOneFloat("maxsampleluminance", Infinity);

        std::string channelnames = params.FindOneString("channels", "R,G,B,A");
        std::string renameas = params.FindOneString("renameas", channelnames);

        return new MultichannelFilm(Point2i(xres, yres), crop, std::move(filter), diagonal,
                        filename, scale, maxSampleLuminance, channelnames, renameas);
    }

}  // namespace pbrt
