#include "scene/light/envmap.h"

#include "core/math.h"
#include "io/image.h"
#include "render/sampling/cdf.h"
#include <cstdio>
#include <tbb/parallel_for.h>
#include <utility>

namespace ybi
{

bool CreateEnvmapCDF(EnvmapLight &light)
{
    uint32_t width = 1;
    uint32_t height = 1;
    std::vector<Vec3f> pixels;
    if (light.texturePath.empty())
    {
        pixels.emplace_back(1.f);
    }
    else
    {
        if (!light.textureFormat.empty() && light.textureFormat != "automatic" &&
            light.textureFormat != "latlong")
        {
            std::fprintf(stderr,
                         "Unsupported environment map projection: %s\n",
                         light.textureFormat.c_str());
            return false;
        }
        if (!LoadLinearRGBImage(light.texturePath, width, height, pixels))
            return false;
        if (light.textureFormat != "latlong" && size_t(width) != 2 * size_t(height))
        {
            std::fprintf(stderr,
                         "Cannot infer environment map projection for %s; "
                         "set the texture format explicitly to latlong.\n",
                         light.texturePath.c_str());
            return false;
        }
    }

    std::vector<float> values(pixels.size());
    tbb::parallel_for(size_t(0), size_t(height), [&](size_t y) {
        const float theta = YBI_PI * (float(y) + 0.5f) / float(height);
        const float sinTheta = Sin(theta);
        for (uint32_t x = 0; x < width; x++)
        {
            const size_t pixelIndex = y * width + x;
            Vec3f &pixel = pixels[pixelIndex];
            pixel *= light.radianceScale;
            const float brightness = (Abs(pixel.x) + Abs(pixel.y) + Abs(pixel.z)) / 3.f;
            values[pixelIndex] = brightness * sinTheta;
        }
    });

    std::vector<float> conditionalCDF;
    std::vector<float> marginalCDF;
    float totalWeight;
    if (!CreatePiecewise2DCDF(values, width, height, conditionalCDF, marginalCDF, totalWeight))
    {
        std::fprintf(
            stderr, "Invalid environment map sampling weights: %s\n", light.texturePath.c_str());
        return false;
    }

    light.width = width;
    light.height = height;
    light.pixels = std::move(pixels);
    light.conditionalCDF = std::move(conditionalCDF);
    light.marginalCDF = std::move(marginalCDF);
    light.totalWeight = totalWeight;
    light.textureFormat = "latlong";
    return true;
}

} // namespace ybi
