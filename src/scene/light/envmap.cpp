#include "scene/light/envmap.h"

#include "core/math.h"
#include "io/image.h"
#include "render/sampling/cdf.h"
#include <atomic>
#include <cmath>
#include <tbb/parallel_for.h>
#include <utility>

namespace ybi
{

static float EnvmapSamplingWeight(const Vec3f &pixel, float sinTheta)
{
    return ((Abs(pixel.x) + Abs(pixel.y) + Abs(pixel.z)) / 3.f) * sinTheta;
}

bool LoadEnvmap(EnvmapLight &light, std::string *error)
{
    if (!light.textureFormat.empty() && light.textureFormat != "automatic" &&
        light.textureFormat != "latlong")
    {
        *error = "Unsupported environment map projection: " + light.textureFormat;
        return false;
    }
    if (!std::isfinite(light.radianceScale.x) || !std::isfinite(light.radianceScale.y) ||
        !std::isfinite(light.radianceScale.z))
    {
        *error = "Non-finite environment map radiance scale: " + light.texturePath;
        return false;
    }

    uint32_t width = 1;
    uint32_t height = 1;
    std::vector<Vec3f> rgbPixels;
    if (light.texturePath.empty())
    {
        rgbPixels.emplace_back(1.f);
    }
    else
    {
        if (!LoadLinearRGBImage(light.texturePath, width, height, rgbPixels))
        {
            *error = "Failed to load environment map: " + light.texturePath;
            return false;
        }
        if (light.textureFormat != "latlong" && size_t(width) != 2 * size_t(height))
        {
            *error = "Cannot infer environment map projection for " + light.texturePath +
                     "; set the texture format explicitly to latlong.";
            return false;
        }
    }

    std::vector<Vec4f> pixels(rgbPixels.size());
    std::vector<float> rowTotals(height);
    std::atomic<bool> valid{true};
    tbb::parallel_for(size_t(0), size_t(height), [&](size_t y) {
        const float theta = YBI_PI * (float(y) + 0.5f) / float(height);
        const float sinTheta = Sin(theta);
        float rowTotal = 0.f;
        for (uint32_t x = 0; x < width; x++)
        {
            const Vec3f pixel = rgbPixels[y * width + x] * light.radianceScale;
            const float weight = EnvmapSamplingWeight(pixel, sinTheta);
            if (!std::isfinite(pixel.x) || !std::isfinite(pixel.y) || !std::isfinite(pixel.z) ||
                !std::isfinite(weight))
            {
                valid.store(false, std::memory_order_relaxed);
                return;
            }
            pixels[y * width + x] = Vec4f(pixel, 1.f);
            rowTotal += weight;
        }
        rowTotals[y] = rowTotal;
    });

    float totalWeight = 0.f;
    for (float rowTotal : rowTotals)
        totalWeight += rowTotal;
    if (!valid.load(std::memory_order_relaxed) || !std::isfinite(totalWeight))
    {
        *error = "Invalid environment map radiance or sampling weights: " + light.texturePath;
        return false;
    }

    light.width = width;
    light.height = height;
    light.pixels = std::move(pixels);
    light.textureFormat = "latlong";
    return true;
}

void PrepareEnvmap(EnvmapLight &light)
{
    const uint32_t width = light.width;
    const uint32_t height = light.height;
    YBI_ASSERT(width > 0 && height > 0);
    YBI_ASSERT(light.pixels.size() == size_t(width) * height);

    std::vector<float> values(light.pixels.size());
    tbb::parallel_for(size_t(0), size_t(height), [&](size_t y) {
        const float theta = YBI_PI * (float(y) + 0.5f) / float(height);
        const float sinTheta = Sin(theta);
        for (uint32_t x = 0; x < width; x++)
        {
            const size_t pixelIndex = y * width + x;
            const Vec3f pixel = light.pixels[pixelIndex].Xyz();
            values[pixelIndex] = EnvmapSamplingWeight(pixel, sinTheta);
        }
    });
    light.totalWeight =
        CreatePiecewise2DCDF(values, width, height, light.conditionalCDF, light.marginalCDF);
}

} // namespace ybi
