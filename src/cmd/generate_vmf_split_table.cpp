#include "core/math.h"
#include "render/path_guiding/lut.h"

#include <array>
#include <cmath>
#include <cstdio>
#include <vector>

namespace ybi
{

struct SolverConfig
{
    double probabilityThreshold;
    double angleThreshold;
    double convergenceThreshold;

    int thetaSteps;
    int phiSteps;
};

double CalculateIntegral(const SolverConfig &config, double invSqrtK, double b)
{
    const double k = invSqrtK == 0.0 ? INFINITY : (1.0 / invSqrtK) / invSqrtK;
    if (std::isinf(k))
    {
        // Infinite-concentration limit.
        return b > config.angleThreshold ? 1.0 : (b < config.angleThreshold ? 0.0 : 0.5);
    }
    if (k == 0.0)
    {
        // Uniform distribution: probability is independent of mean direction.
        return 0.5 * (1.0 + Cos(config.angleThreshold));
    }

    const double dphi = 2.0 * YBI_PI / config.phiSteps;
    const double dtheta = (YBI_PI - config.angleThreshold) / config.thetaSteps;
    const double meanX = Sin(b);
    const double meanZ = Cos(b);
    double integral = 0.0;
    std::vector<double> cosPhi(config.phiSteps);
    for (int phiIndex = 0; phiIndex < config.phiSteps; ++phiIndex)
    {
        cosPhi[phiIndex] = Cos((phiIndex + 0.5) * dphi);
    }

    for (int thetaIndex = 0; thetaIndex < config.thetaSteps; ++thetaIndex)
    {
        const double theta = config.angleThreshold + (thetaIndex + 0.5) * dtheta;
        const double sinTheta = Sin(theta);
        const double cosTheta = Cos(theta);
        double azimuthIntegral = 0.0;
        for (int phiIndex = 0; phiIndex < config.phiSteps; ++phiIndex)
        {
            const double cosine = sinTheta * cosPhi[phiIndex] * meanX + cosTheta * meanZ;
            // Prevent rounding from creating a positive exponent at large k.
            const double exponent = k * (std::fmin(1.0, cosine) - 1.0);
            azimuthIntegral += std::exp(exponent);
        }
        integral += azimuthIntegral * sinTheta;
    }

    const double normalization = k / (2.0 * YBI_PI * -std::expm1(-2.0 * k));
    return integral * dtheta * dphi * normalization;
}

// One bisection step; true means converged. The caller supplies a valid bracket.
bool Bisection(const SolverConfig &config,
               double angle,
               double angleLow,
               double angleHigh,
               double value,
               double valueLow,
               double &newAngle,
               double &newAngleLow,
               double &newAngleHigh)
{
    const double f = value - config.probabilityThreshold;
    const double fLow = valueLow - config.probabilityThreshold;
    newAngle = angle;
    newAngleLow = angleLow;
    newAngleHigh = angleHigh;

    if (std::abs(f) <= config.convergenceThreshold)
    {
        return true;
    }

    if ((f < 0.0) != (fLow < 0.0))
    {
        newAngleHigh = angle;
    }
    else
    {
        newAngleLow = angle;
    }
    newAngle = newAngleLow + 0.5 * (newAngleHigh - newAngleLow);
    // Stop when no representable midpoint remains.
    if (newAngle == newAngleLow || newAngle == newAngleHigh)
    {
        newAngle = newAngleHigh;
        return true;
    }
    return false;
}

// Closed-form probabilities for beta=0 and beta=pi.
double CalculateEndpointProbability(const SolverConfig &config, double invSqrtK, double beta)
{
    if (invSqrtK == 0.0)
    {
        return beta == 0.0 ? 0.0 : 1.0;
    }
    const double k = (1.0 / invSqrtK) / invSqrtK;
    const double c = Cos(config.angleThreshold);
    const double denominator = -std::expm1(-2.0 * k);
    const double capFactor = -std::expm1(-k * (1.0 + c));
    return beta == 0.0 ? std::exp(k * (c - 1.0)) * capFactor / denominator
                       : capFactor / denominator;
}

double FindMaxInvSqrtK(const SolverConfig &config, double &endpointBeta)
{
    const double uniformProbability = 0.5 * (1.0 + Cos(config.angleThreshold));
    endpointBeta = uniformProbability < config.probabilityThreshold ? YBI_PI : 0.0;
    double low = 0.0;
    double high = 1.0;
    const double valueLow = CalculateEndpointProbability(config, low, endpointBeta);
    while ((CalculateEndpointProbability(config, high, endpointBeta) <
            config.probabilityThreshold) == (valueLow < config.probabilityThreshold))
    {
        high *= 2.0;
    }

    double invSqrtK = 0.5 * (low + high);
    for (;;)
    {
        double next;
        const bool done = Bisection(config,
                                    invSqrtK,
                                    low,
                                    high,
                                    CalculateEndpointProbability(config, invSqrtK, endpointBeta),
                                    CalculateEndpointProbability(config, low, endpointBeta),
                                    next,
                                    low,
                                    high);
        if (done)
        {
            return next;
        }
        invSqrtK = next;
    }
}

double FindSplitAngle(SolverConfig config, double invSqrtK)
{
    if (invSqrtK == 0.0)
    {
        return config.angleThreshold;
    }

    const double k = (1.0 / invSqrtK) / invSqrtK;
    const double alpha = 1.0 - config.probabilityThreshold;
    // The 1-alpha confidence cone supplies an upper bracket for beta.
    const double coneCosine =
        1.0 + std::log(std::exp(-2.0 * k) + alpha * -std::expm1(-2.0 * k)) / k;
    double low = 0.0;
    double high =
        Min(double(YBI_PI), config.angleThreshold + std::acos(Clamp(coneCosine, -1.0, 1.0)));
    config.thetaSteps =
        Max(config.thetaSteps, int(std::ceil(32.0 * (YBI_PI - config.angleThreshold) / invSqrtK)));
    config.phiSteps =
        Max(config.phiSteps,
            int(std::ceil(4.0 * YBI_PI * Sin(Min(high, double(YBI_PI) / 2.0)) / invSqrtK)));

    double valueLow = CalculateEndpointProbability(config, invSqrtK, low);
    const double valueHigh = high == YBI_PI ? CalculateEndpointProbability(config, invSqrtK, high)
                                            : CalculateIntegral(config, invSqrtK, high);
    if (valueLow >= config.probabilityThreshold)
    {
        return 0.0;
    }
    if (valueHigh == config.probabilityThreshold)
    {
        return high;
    }

    double angle = 0.5 * (low + high);
    for (;;)
    {
        const double value = CalculateIntegral(config, invSqrtK, angle);
        double next, nextLow, nextHigh;
        if (Bisection(config, angle, low, high, value, valueLow, next, nextLow, nextHigh))
        {
            return next;
        }
        if (nextLow != low)
        {
            valueLow = value;
        }
        low = nextLow;
        high = nextHigh;
        angle = next;
    }
}

} // namespace ybi

int main(int argc, char **argv)
{
    using namespace ybi;
    if (argc > 2)
    {
        std::fprintf(stderr, "Usage: %s [output-file]\n", argv[0]);
        return 1;
    }
    const char *outputPath = argc == 2 ? argv[1] : VMFSplitTableFilename;
    FILE *output = std::fopen(outputPath, "wb");
    if (!output)
    {
        std::fprintf(stderr, "Failed to open output file: %s\n", outputPath);
        return 1;
    }
    constexpr int numEntries = VMFSplitTableSize;
    // Paper defaults: alpha=1e-4 and tau=3 degrees.
    const SolverConfig config{1.0 - 1e-4, 3.0 * YBI_PI / 180.0, 1e-10, 4096, 64};
    double endpointBeta;
    const double maxInvSqrtK = FindMaxInvSqrtK(config, endpointBeta);
    std::array<float, numEntries + 1> values;
    values[0] = static_cast<float>(maxInvSqrtK);
    std::printf("k,inv_sqrt_k,cos_beta\n");
    for (int index = 0; index < numEntries; ++index)
    {
        const double invSqrtK = maxInvSqrtK * index / (numEntries - 1);
        const double beta =
            index == numEntries - 1 ? endpointBeta : FindSplitAngle(config, invSqrtK);
        const double cosBeta = index == numEntries - 1 ?
            (endpointBeta == 0.0 ? 1.0 : -1.0) : Cos(beta);
        values[index + 1] = static_cast<float>(cosBeta);
        std::printf("%.17g, %.17g,%.17g\n", 1 / Sqr(invSqrtK), invSqrtK, cosBeta);
    }
    const size_t written = std::fwrite(values.data(), sizeof(float), values.size(), output);
    const int closeResult = std::fclose(output);
    if (written != values.size() || closeResult != 0)
    {
        std::fprintf(stderr, "Failed to write output file: %s\n", outputPath);
        return 1;
    }
    std::fprintf(stderr, "Wrote %u entries (%zu bytes) to %s\n",
                 VMFSplitTableSize, values.size() * sizeof(float), outputPath);
    return 0;
}
