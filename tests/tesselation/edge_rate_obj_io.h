#pragma once

#include <pxr/base/gf/vec3f.h>

#include <string>
#include <vector>

namespace ybi::testio
{
struct EdgeRateDebugLine
{
    std::vector<pxr::GfVec3f> points;
    int rate = 1;
};

bool WriteEdgeRateDebugObj(const std::string &path,
                           const std::vector<EdgeRateDebugLine> &lines,
                           int minRate,
                           int maxRate);
} // namespace ybi::testio

