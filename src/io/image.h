#pragma once

#include "core/vec3f.h"
#include "device/device_texture.h"
#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

namespace ybi
{

struct Image
{
    void *ptr = nullptr;
    int width = 0;
    int height = 0;
    int numChannels = 0;
    int numBytesPerChannel = 0;
};

int GetImageFileChannelCount(const std::string &filename);
SampleFormat LoadImageFile(Image *image, const std::string &filename);
void FlipImageVertically(Image *image);
float ReadChannelAsFloat(const void *data,
                         size_t pixelIndex,
                         int channel,
                         int numChannels,
                         int numBytesPerChannel,
                         SampleFormat sampleFormat);
void WriteChannelFromFloat(void *data,
                           size_t pixelIndex,
                           int channel,
                           int numChannels,
                           int numBytesPerChannel,
                           SampleFormat sampleFormat,
                           float value);
void PadRGBToRGBA(Image *image, SampleFormat sampleFormat);
void LinearizeSrgbTexture(Image *image, SampleFormat sampleFormat);

bool LoadLinearRGBImage(const std::string &filename,
                        uint32_t &width,
                        uint32_t &height,
                        std::vector<Vec3f> &pixels);

} // namespace ybi
