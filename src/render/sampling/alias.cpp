#include "render/sampling/alias.h"

#include "core/assert.h"
#include <tbb/parallel_for.h>

namespace ybi
{

float CreateAliasTable(const float *weights, uint32_t count, AliasEntry *table)
{
    YBI_ASSERT(count != 0);

    float total = 0.f;
    for (uint32_t i = 0; i < count; i++)
        total += weights[i];

    std::vector<uint32_t> small;
    std::vector<uint32_t> large;
    small.reserve(count);
    large.reserve(count);
    for (uint32_t i = 0; i < count; i++)
    {
        const float probability = total == 0.f ? 1.f / float(count) : weights[i] / total;
        table[i].probability = probability;
        table[i].threshold = probability * float(count);
        table[i].alias = i;
        if (table[i].threshold < 1.f)
            small.push_back(i);
        else
            large.push_back(i);
    }

    while (!small.empty() && !large.empty())
    {
        const uint32_t s = small.back();
        small.pop_back();
        const uint32_t l = large.back();
        large.pop_back();

        table[s].alias = l;
        table[l].threshold = (table[l].threshold - 1.f) + table[s].threshold;
        if (table[l].threshold < 1.f)
            small.push_back(l);
        else
            large.push_back(l);
    }

    for (uint32_t i : small)
        table[i].threshold = 1.f;
    for (uint32_t i : large)
        table[i].threshold = 1.f;
    return total;
}

float CreateAliasTable2D(const std::vector<float> &weights,
                         uint32_t width,
                         uint32_t height,
                         std::vector<AliasEntry> &conditionalAlias,
                         std::vector<AliasEntry> &marginalAlias)
{
    YBI_ASSERT(width > 0 && height > 0);
    YBI_ASSERT(weights.size() == size_t(width) * height);

    conditionalAlias.resize(size_t(width) * height);
    marginalAlias.resize(height);
    std::vector<float> rowTotals(height);
    tbb::parallel_for(size_t(0), size_t(height), [&](size_t y) {
        rowTotals[y] = CreateAliasTable(weights.data() + y * width,
                                        width,
                                        conditionalAlias.data() + y * width);
    });
    return CreateAliasTable(rowTotals.data(), height, marginalAlias.data());
}

} // namespace ybi
