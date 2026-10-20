#include "scene/light/envmap.h"

#include "core/math.h"
#include "io/image.h"
#include "render/sampling/alias.h"
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

void PrepareEnvmap(EnvmapLight &light,
                    uint32_t maxDistributionWidth,
                    uint32_t maxDistributionHeight)
{
    const uint32_t width = light.width;
    const uint32_t height = light.height;
    YBI_ASSERT(width > 0 && height > 0);
    YBI_ASSERT(light.pixels.size() == size_t(width) * height);
    YBI_ASSERT(maxDistributionWidth > 0 && maxDistributionHeight > 0);

    const uint32_t distributionWidth = Min(width, maxDistributionWidth);
    const uint32_t distributionHeight = Min(height, maxDistributionHeight);
    std::vector<float> values(size_t(distributionWidth) * distributionHeight);
    std::vector<float> sinTheta(height);
    for (uint32_t y = 0; y < height; y++)
        sinTheta[y] = Sin(YBI_PI * (float(y) + 0.5f) / float(height));

    tbb::parallel_for(size_t(0), size_t(distributionHeight), [&](size_t y) {
        const float y0 = float(y) * float(height) / float(distributionHeight);
        const float y1 = float(y + 1) * float(height) / float(distributionHeight);
        const uint32_t sourceYBegin = static_cast<uint32_t>(Floor(y0));
        const uint32_t sourceYEnd = Min(static_cast<uint32_t>(ceilf(y1)), height);
        for (uint32_t x = 0; x < distributionWidth; x++)
        {
            const float x0 = float(x) * float(width) / float(distributionWidth);
            const float x1 = float(x + 1) * float(width) / float(distributionWidth);
            const uint32_t sourceXBegin = static_cast<uint32_t>(Floor(x0));
            const uint32_t sourceXEnd = Min(static_cast<uint32_t>(ceilf(x1)), width);
            float weight = 0.f;
            for (uint32_t sourceY = sourceYBegin; sourceY < sourceYEnd; sourceY++)
            {
                const float overlapY =
                    Min(y1, float(sourceY + 1)) - Max(y0, float(sourceY));
                float rowWeight = 0.f;
                for (uint32_t sourceX = sourceXBegin; sourceX < sourceXEnd; sourceX++)
                {
                    const float overlapX =
                        Min(x1, float(sourceX + 1)) - Max(x0, float(sourceX));
                    const Vec3f pixel = light.pixels[size_t(sourceY) * width + sourceX].Xyz();
                    rowWeight += overlapX * EnvmapSamplingWeight(pixel, sinTheta[sourceY]);
                }
                weight += overlapY * rowWeight;
            }
            values[y * distributionWidth + x] = weight / ((x1 - x0) * (y1 - y0));
        }
    });
    light.totalWeight = CreateAliasTable2D(values, distributionWidth, distributionHeight,
                                           light.conditionalAlias, light.marginalAlias);
    light.distributionWidth = distributionWidth;
    light.distributionHeight = distributionHeight;
}

} // namespace ybi
