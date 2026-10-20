#pragma once

#include "core/math.h"
#include "core/vec2f.h"
#include "kernel/sampling/types.h"

namespace ybi
{

YBI_CPU_GPU uint32_t SampleAliasTable(
    const AliasEntry *table, uint32_t count, float u, float &fraction, float &probability)
{
    const float position = Min(u, FLOAT_BEFORE_ONE) * float(count);
    const uint32_t bin = Min(static_cast<uint32_t>(position), count - 1);
    const float offset = Min(position - float(bin), FLOAT_BEFORE_ONE);
    const AliasEntry &entry = table[bin];

    uint32_t index;
    if (offset < entry.threshold)
    {
        index = bin;
        fraction = offset / entry.threshold;
    }
    else
    {
        index = entry.alias;
        fraction = (offset - entry.threshold) / (1.f - entry.threshold);
    }
    fraction = Min(fraction, FLOAT_BEFORE_ONE);
    probability = table[index].probability;
    YBI_ASSERT(probability != 0.f);
    return index;
}

YBI_CPU_GPU bool SampleAliasTable2D(const AliasEntry *conditionalAlias,
                                    const AliasEntry *marginalAlias,
                                    uint32_t width,
                                    uint32_t height,
                                    const Vec2f &u,
                                    Vec2f &uv,
                                    float &pdf)
{
    float dx, dy;
    float rowProbability, columnProbability;
    const uint32_t y = SampleAliasTable(marginalAlias, height, u.y, dy, rowProbability);
    const AliasEntry *row = conditionalAlias + size_t(y) * width;
    const uint32_t x = SampleAliasTable(row, width, u.x, dx, columnProbability);

    uv = Vec2f(Min((float(x) + dx) / float(width), FLOAT_BEFORE_ONE),
               Min((float(y) + dy) / float(height), FLOAT_BEFORE_ONE));
    pdf = rowProbability * columnProbability * float(width) * float(height);
    return pdf > 0.f;
}

YBI_CPU_GPU float PDFAliasTable2D(const AliasEntry *conditionalAlias,
                                  const AliasEntry *marginalAlias,
                                  uint32_t width,
                                  uint32_t height,
                                  const Vec2f &uv)
{
    const uint32_t x = Min(static_cast<uint32_t>(uv.x * width), width - 1);
    const uint32_t y = Min(static_cast<uint32_t>(uv.y * height), height - 1);
    return marginalAlias[y].probability * conditionalAlias[size_t(y) * width + x].probability *
           float(width) * float(height);
}

} // namespace ybi
