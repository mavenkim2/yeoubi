#pragma once

#include "device/cuda_compat.h"

namespace ybi
{

YBI_CPU_GPU float PowerHeuristic(float sampledPdf, float otherPdf)
{
    const float sampledPdfSquared = sampledPdf * sampledPdf;
    const float otherPdfSquared = otherPdf * otherPdf;
    const float sum = sampledPdfSquared + otherPdfSquared;
    return sum > 0.f ? sampledPdfSquared / sum : 0.f;
}

} // namespace ybi
