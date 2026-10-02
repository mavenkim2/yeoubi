#include "kernel/bsdf/openpbr/ggx.h"
#include "kernel/integrator/rng.h"
#include "kernel/integrator/sampling.h"

#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <vector>

using namespace ybi;

namespace
{

constexpr int TableSize = 32;
constexpr int DefaultSampleCount = 10000;
constexpr const char *DefaultOutputPath = "openpbr_opaque_dielectric_energy_complement.bin";
constexpr const char *DefaultAverageOutputPath =
    "openpbr_opaque_dielectric_average_energy_complement.bin";
constexpr float IorMax = 2.5f;
constexpr float MinIorAboveOne = 1.0001f;

enum class TableMode
{
    Directional,
    Albedo,
    All,
};

// Match OpenPBR's LUT parameterization. Half of the IOR axis covers internal
// reflection from 1 / IorMax to 1, and the other half covers 1 to IorMax.
float IorFromTableIndex(int index)
{
    constexpr int halfTableSize = TableSize / 2;
    constexpr int halfTableSizeMinusOne = halfTableSize - 1;

    // The two slices adjacent to IOR 1 use the same safety clamp as OpenPBR's
    // refraction setup. This avoids the singular indexed-matched grazing case.
    if (index == halfTableSizeMinusOne)
        return 1.f / MinIorAboveOne;
    if (index == halfTableSize)
        return MinIorAboveOne;

    if (index < halfTableSize)
    {
        const float fraction = float(halfTableSizeMinusOne - index) / float(halfTableSizeMinusOne);
        return 1.f / (1.f + fraction * (IorMax - 1.f));
    }

    const float fraction = float(index - halfTableSize) / float(halfTableSizeMinusOne);
    return 1.f + fraction * (IorMax - 1.f);
}

float EstimateDirectionalReflectance(
    float ior, float alpha, float cosTheta, int sampleCount, uint32_t seed)
{
    // Indexed-matched media have no interface reflection, including at grazing.
    if (ior == 1.f)
        return 0.f;

    // SampleGGX requires a direction strictly above the macrosurface. The zero
    // cosine table entry stores the grazing-angle limit.
    const float sampleCosTheta = Max(cosTheta, 1.0e-3f);

    // GGX becomes a delta distribution at alpha zero, so use its analytic limit.
    if (alpha == 0.f)
        return FrDielectric(sampleCosTheta, ior);

    const float sinTheta = sqrtf(Max(0.f, 1.f - sampleCosTheta * sampleCosTheta));
    const Vec3f viewDirection(sinTheta, 0.f, sampleCosTheta);
    const Vec2f roughness(alpha);

    RNG rng = RNG::Init(seed);
    double reflectanceSum = 0.0;
    for (int sampleIndex = 0; sampleIndex < sampleCount; ++sampleIndex)
    {
        Vec3f lightDirection;
        Vec3f weight;
        if (SampleGGX(viewDirection, roughness, rng.Uniform2D(), ior, lightDirection, weight))
        {
            // VNDF sampling makes weight equal to f_r * cos(theta_l) / pdf.
            reflectanceSum += double(weight.x);
        }
    }

    return float(reflectanceSum / double(sampleCount));
}

float EstimateAverageReflectance(float ior, float alpha, int sampleCount, uint32_t seed)
{
    if (ior == 1.f)
        return 0.f;

    const Vec2f roughness(alpha);
    RNG rng = RNG::Init(seed);
    double reflectanceSum = 0.0;
    for (int sampleIndex = 0; sampleIndex < sampleCount; ++sampleIndex)
    {
        // The hemispherical-hemispherical reflectance is
        //   (1 / pi) integral integral f(wi, wo) cos_i cos_o dwi dwo.
        // Sampling wo with p(wo) = cos_o / pi cancels both the outer
        // cosine and the 1 / pi normalization. SampleGGX's returned weight
        // is f(wi, wo) * cos_i / p(wi | wo), so it is the complete estimator.
        const Vec3f viewDirection = SampleCosineHemisphere(rng.Uniform2D());

        if (alpha == 0.f)
        {
            reflectanceSum += double(FrDielectric(viewDirection.z, ior));
            continue;
        }

        Vec3f lightDirection;
        Vec3f weight;
        if (SampleGGX(viewDirection, roughness, rng.Uniform2D(), ior, lightDirection, weight))
            reflectanceSum += double(weight.x);
    }

    return float(reflectanceSum / double(sampleCount));
}

bool WriteValues(const char *outputPath, const std::vector<uint16_t> &values)
{
    FILE *output = std::fopen(outputPath, "wb");
    if (!output)
    {
        std::fprintf(stderr, "Failed to open output file: %s\n", outputPath);
        return false;
    }

    const size_t written = std::fwrite(values.data(), sizeof(uint16_t), values.size(), output);
    const int closeResult = std::fclose(output);
    if (written != values.size() || closeResult != 0)
    {
        std::fprintf(stderr, "Failed to write output file: %s\n", outputPath);
        return false;
    }

    std::printf("Wrote %zu uint16 values (%zu bytes) to %s\n",
                values.size(),
                values.size() * sizeof(uint16_t),
                outputPath);
    return true;
}

int ParseSampleCount(const char *text)
{
    char *end = nullptr;
    const long value = std::strtol(text, &end, 10);
    if (end == text || *end != '\0' || value <= 0)
        return 0;

    return int(value);
}

bool ParseTableMode(const char *text, TableMode &mode)
{
    if (std::strcmp(text, "directional") == 0)
        mode = TableMode::Directional;
    else if (std::strcmp(text, "albedo") == 0)
        mode = TableMode::Albedo;
    else if (std::strcmp(text, "all") == 0)
        mode = TableMode::All;
    else
        return false;

    return true;
}

void PrintUsage(const char *executablePath)
{
    std::fprintf(stderr,
                 "Usage:\n"
                 "  %s directional [samples-per-entry] [output-file]\n"
                 "  %s albedo     [samples-per-entry] [output-file]\n"
                 "  %s all        [samples-per-entry] [directional-output-file] "
                 "[albedo-output-file]\n",
                 executablePath,
                 executablePath,
                 executablePath);
}

} // namespace

int main(int argc, char **argv)
{
    TableMode mode;
    if (argc < 2 || !ParseTableMode(argv[1], mode))
    {
        PrintUsage(argv[0]);
        return 1;
    }

    const bool generateDirectional = mode == TableMode::Directional || mode == TableMode::All;
    const bool generateAlbedo = mode == TableMode::Albedo || mode == TableMode::All;
    const int maximumArgumentCount = mode == TableMode::All ? 5 : 4;
    if (argc > maximumArgumentCount)
    {
        PrintUsage(argv[0]);
        return 1;
    }

    const int sampleCount = argc >= 3 ? ParseSampleCount(argv[2]) : DefaultSampleCount;
    if (sampleCount == 0)
    {
        std::fprintf(stderr, "samples-per-entry must be a positive integer.\n");
        return 1;
    }

    const char *directionalOutputPath = DefaultOutputPath;
    const char *albedoOutputPath = DefaultAverageOutputPath;
    if (mode == TableMode::Directional && argc >= 4)
        directionalOutputPath = argv[3];
    else if (mode == TableMode::Albedo && argc >= 4)
        albedoOutputPath = argv[3];
    else if (mode == TableMode::All)
    {
        if (argc >= 4)
            directionalOutputPath = argv[3];
        if (argc >= 5)
            albedoOutputPath = argv[4];
    }

    if (generateDirectional)
    {
        std::vector<uint16_t> values;
        values.reserve(TableSize * TableSize * TableSize);

        uint32_t cellIndex = 0;
        for (int iorIndex = 0; iorIndex < TableSize; ++iorIndex)
        {
            const float ior = IorFromTableIndex(iorIndex);
            for (int alphaIndex = 0; alphaIndex < TableSize; ++alphaIndex)
            {
                // The LUT is uniform in perceptual roughness and stores alpha = roughness^2.
                const float perceptualRoughness = float(alphaIndex) / float(TableSize - 1);
                const float alpha = perceptualRoughness * perceptualRoughness;

                for (int cosThetaIndex = 0; cosThetaIndex < TableSize; ++cosThetaIndex)
                {
                    const float cosTheta = float(cosThetaIndex) / float(TableSize - 1);
                    const float directionalReflectance = EstimateDirectionalReflectance(
                        ior, alpha, cosTheta, sampleCount, 0x9e3779b9u + cellIndex);
                    const float energyComplement =
                        Clamp(1.f - directionalReflectance, 0.f, 1.f);

                    values.push_back(
                        static_cast<uint16_t>(std::lround(energyComplement * 65535.f)));
                    ++cellIndex;
                }
            }
        }

        if (!WriteValues(directionalOutputPath, values))
            return 1;
    }

    if (generateAlbedo)
    {
        std::vector<uint16_t> averageValues;
        averageValues.reserve(TableSize * TableSize);

        uint32_t averageCellIndex = 0;
        for (int iorIndex = 0; iorIndex < TableSize; ++iorIndex)
        {
            const float ior = IorFromTableIndex(iorIndex);
            for (int alphaIndex = 0; alphaIndex < TableSize; ++alphaIndex)
            {
                const float perceptualRoughness = float(alphaIndex) / float(TableSize - 1);
                const float alpha = perceptualRoughness * perceptualRoughness;
                const float averageReflectance = EstimateAverageReflectance(
                    ior, alpha, sampleCount, 0x85ebca6bu + averageCellIndex);
                const float averageEnergyComplement =
                    Clamp(1.f - averageReflectance, 0.f, 1.f);

                averageValues.push_back(
                    static_cast<uint16_t>(std::lround(averageEnergyComplement * 65535.f)));
                ++averageCellIndex;
            }
        }

        if (!WriteValues(albedoOutputPath, averageValues))
            return 1;
    }

    return 0;
}
