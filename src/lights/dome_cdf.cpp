#include "lights/dome_cdf.h"

#include "util/cdf.h"
#include "util/half_float.h"
#include "util/math_constants.h"
#include "util/math_common.h"
#include "util/vec3.h"

#include <algorithm>
#include <cmath>
#include <cstring>

namespace ybi
{
namespace lights
{
namespace
{

Vec3 ReadDecodedTextureRgb(const DecodedMaterialTexture &texture, int x, int y)
{
    if (!texture.valid || texture.width <= 0 || texture.height <= 0 || texture.pixels.empty())
    {
        return Vec3(0.0f);
    }

    x = Clamp(x, 0, texture.width - 1);
    y = Clamp(y, 0, texture.height - 1);
    const size_t pixelIndex =
        static_cast<size_t>(y) * static_cast<size_t>(texture.width) + static_cast<size_t>(x);

    switch (texture.format)
    {
        case TextureFormat::RGBA8_UNORM:
        {
            const size_t offset = pixelIndex * 4u;
            if (offset + 3u >= texture.pixels.size())
            {
                return Vec3(0.0f);
            }
            return Vec3(float(texture.pixels[offset + 0u]) * (1.0f / 255.0f),
                        float(texture.pixels[offset + 1u]) * (1.0f / 255.0f),
                        float(texture.pixels[offset + 2u]) * (1.0f / 255.0f));
        }
        case TextureFormat::RGBA16_FLOAT:
        {
            const size_t byteOffset = pixelIndex * 4u * sizeof(uint16_t);
            if (byteOffset + 4u * sizeof(uint16_t) > texture.pixels.size())
            {
                return Vec3(0.0f);
            }

            uint16_t rBits = 0u;
            uint16_t gBits = 0u;
            uint16_t bBits = 0u;
            std::memcpy(&rBits,
                        texture.pixels.data() + byteOffset + 0u * sizeof(uint16_t),
                        sizeof(uint16_t));
            std::memcpy(&gBits,
                        texture.pixels.data() + byteOffset + 1u * sizeof(uint16_t),
                        sizeof(uint16_t));
            std::memcpy(&bBits,
                        texture.pixels.data() + byteOffset + 2u * sizeof(uint16_t),
                        sizeof(uint16_t));
            return Vec3(util::HalfBitsToFloat(rBits),
                        util::HalfBitsToFloat(gBits),
                        util::HalfBitsToFloat(bBits));
        }
        default:
            return Vec3(0.0f);
    }
}

} // namespace

bool BuildDomeTextureCdf(const DecodedMaterialTexture &domeTexture,
                         std::vector<float> *outConditional,
                         std::vector<float> *outMarginal,
                         std::string *outError)
{
    YBI_ASSERT(outConditional);
    YBI_ASSERT(outMarginal);
    outConditional->clear();
    outMarginal->clear();

    if (!domeTexture.valid || domeTexture.width <= 0 || domeTexture.height <= 0)
    {
        return true;
    }

    if (domeTexture.format != TextureFormat::RGBA8_UNORM &&
        domeTexture.format != TextureFormat::RGBA16_FLOAT)
    {
        if (outError)
        {
            *outError = "unsupported dome texture format for CDF";
        }
        return false;
    }

    outConditional->resize(static_cast<size_t>(domeTexture.height) *
                           static_cast<size_t>(domeTexture.width + 1));
    outMarginal->resize(static_cast<size_t>(domeTexture.height + 1));

    InitializeCDF2D(
        [&](float u, float v) {
            const int x = Clamp(int(u * float(domeTexture.width)), 0, domeTexture.width - 1);
            const int y = Clamp(int(v * float(domeTexture.height)), 0, domeTexture.height - 1);
            const Vec3 rgb = ReadDecodedTextureRgb(domeTexture, x, y);
            const float theta = kPi * ((float(y) + 0.5f) / float(domeTexture.height));
            return MaxComponent(rgb) > 0.0f ? Luminance(rgb) * std::max(sinf(theta), 0.0f) : 0.0f;
        },
        Vec2(0.0f),
        Vec2(1.0f),
        domeTexture.width,
        domeTexture.height,
        outConditional->data(),
        outMarginal->data());

    return true;
}

} // namespace lights
} // namespace ybi
