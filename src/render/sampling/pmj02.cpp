#include "render/sampling/pmj02.h"
#include "core/array.h"
#include "core/assert.h"
#include "core/bit.h"
#include "core/vec2f.h"
#include "core/vec2u.h"
#include <vector>

namespace ybi
{

// For a given class of elementary intervals (e.g 1x16 or 2x8), calculate 2D index
// where it belongs.
// e.g class 2x8, (0.06, 0.13) -> (0, 1)
static Vec2u CalculateElementaryIntervalIndex(const Vec2f &sample,
                                              uint32_t intervalIndex,
                                              uint32_t log2SamplesInIteration)
{
    YBI_ASSERT(intervalIndex < log2SamplesInIteration);

    const uint32_t xSlices = 1u << intervalIndex;
    const uint32_t ySlices = 1u << log2SamplesInIteration - intervalIndex;
    const Vec2f interval = Floor(sample * Vec2f(float(xSlices), float(ySlices)));

    YBI_ASSERT(uint32_t(interval.x) < xSlices);
    YBI_ASSERT(uint32_t(interval.y) < ySlices);

    return Vec2u(static_cast<uint32_t>(interval.x), static_cast<uint32_t>(interval.y));
}

static uint32_t Get1DIndex(const Vec2u &index, uint32_t intervalIndex)
{
    return index.y * (1u << intervalIndex) + index.x;
}

static uint32_t CalculateIntervalIndex1D(const Vec2f &sample,
                                         uint32_t intervalIndex,
                                         uint32_t log2SamplesInIteration)
{
    const Vec2u index =
        CalculateElementaryIntervalIndex(sample, intervalIndex, log2SamplesInIteration);
    return Get1DIndex(index, intervalIndex);
}

// Implementation of Efficient Generation of Points that Satisfy
// Two-Dimensional Elementary Intervals
#if 0
static void GetValidXOffsets(int x,
                             int y,
                             int intervalIndex,
                             const std::vector<std::vector<bool>> &occupancy,
                             std::vector<int> &xOffsets,
                             std::vector<int> &yOffsets,
                             int step = 0)
{
    const int xSlices = 1u << intervalIndex;
    const bool occupied = occupancy[intervalIndex][y * xSlices + x];

    if (!occupied)
    {
        if (intervalIndex == static_cast<int>(occupancy.size() - 1))
        {
            xOffsets.push_back(x);
        }
        else
        {
            GetValidXOffsets(x * 2, y >> 1, intervalIndex + 1, occupancy, xOffsets);
            GetValidXOffsets(x * 2 + 1, y >> 1, intervalIndex + 1, occupancy, xOffsets);
        }
    }
}

static void GetValidYOffsets(int x,
                             int y,
                             int intervalIndex,
                             const std::vector<std::vector<bool>> &occupancy,
                             std::vector<int> &yOffsets)
{
    const int xSlices = 1u << intervalIndex;
    const bool occupied = occupancy[intervalIndex][y * xSlices + x];

    if (!occupied)
    {
        if (intervalIndex == 0)
        {
            xOffsets.push_back(x);
        }
        else
        {
            GetValidYOffsets(x >> 1, y * 2 + 1, intervalIndex - 1, occupancy, yOffsets);
            GetValidYOffsets(x >> 1, y * 2 + 1, intervalIndex - 1, occupancy, yOffsets);
        }
    }
}
#endif

void Something()
{
    uint32_t n = 0;
    uint32_t numSamples = 0;

    std::vector<Vec2f> samples;

    int level = 1;
    while (n < numSamples)
    {
        std::vector<std::vector<bool>> elementaryIntervalOccupancy;
        const uint32_t samplesInIteration = 1u << (level * 2);
        const uint32_t log2SamplesInIteration = FindMostSignificantBit(samplesInIteration);

        elementaryIntervalOccupancy.resize(log2SamplesInIteration + 1);
        for (std::vector<bool> &bitVector : elementaryIntervalOccupancy)
        {
            bitVector.assign(samplesInIteration, false);
        }

        // Set the global occupancy
        for (const Vec2f &sample : samples)
        {
            for (size_t index = 0; index < elementaryIntervalOccupancy.size(); index++)
            {
                std::vector<bool> &bitVector = elementaryIntervalOccupancy[index];
                const uint32_t intervalIndex =
                    CalculateIntervalIndex1D(sample, index, log2SamplesInIteration);
                bitVector[intervalIndex] = true;
            }
        }

        // Generate samples in the opposite sub quadrant
        const size_t prevLevelEnd = samples.size();

        // index into elementaryIntervalOccupancy for nxn elementary intervals
        const uint32_t squareIndex = log2SamplesInIteration >> 1;

        for (size_t i = 0; i < prevLevelEnd; i++)
        {
            const Vec2f &sample = samples[i];

            Vec2u intervalIndex =
                CalculateElementaryIntervalIndex(sample, squareIndex, log2SamplesInIteration);

            // Flip to opposite quadrant
            intervalIndex ^= 1;

            const uint32_t oppositeQuadrantIndex = Get1DIndex(intervalIndex, squareIndex);
            YBI_ASSERT(oppositeQuadrantIndex < samplesInIteration);
            YBI_ASSERT(elementaryIntervalOccupancy[squareIndex][oppositeQuadrantIndex] == false);
        }

        // Randomly select one of the two remaining quadrants
        for (size_t i = 0; i < prevLevelEnd; i++)
        {
        }

        level++;
    }
}
} // namespace ybi
