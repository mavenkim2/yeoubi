#include "tesselation/edge_rate_obj_io.h"

#include <algorithm>
#include <fstream>

namespace ybi::testio
{
namespace
{
pxr::GfVec3f RateToColor(int rate, int minRate, int maxRate)
{
    if (maxRate <= minRate)
    {
        return pxr::GfVec3f(0.0f, 1.0f, 0.0f);
    }
    const float t = std::clamp(float(rate - minRate) / float(maxRate - minRate), 0.0f, 1.0f);
    return pxr::GfVec3f(t, 1.0f - t, 0.0f);
}

float ComputeDebugWidth(const std::vector<EdgeRateDebugLine> &lines)
{
    pxr::GfVec3f bbMin(1e30f, 1e30f, 1e30f);
    pxr::GfVec3f bbMax(-1e30f, -1e30f, -1e30f);
    bool any = false;
    for (const EdgeRateDebugLine &line : lines)
    {
        for (const pxr::GfVec3f &p : line.points)
        {
            bbMin[0] = std::min(bbMin[0], p[0]);
            bbMin[1] = std::min(bbMin[1], p[1]);
            bbMin[2] = std::min(bbMin[2], p[2]);
            bbMax[0] = std::max(bbMax[0], p[0]);
            bbMax[1] = std::max(bbMax[1], p[1]);
            bbMax[2] = std::max(bbMax[2], p[2]);
            any = true;
        }
    }
    if (!any)
    {
        return 1e-3f;
    }
    const float diag = (bbMax - bbMin).GetLength();
    return std::max(diag * 0.0025f, 2.5e-4f);
}

pxr::GfVec3f SafePerp(const pxr::GfVec3f &tangent)
{
    pxr::GfVec3f axis(0.0f, 1.0f, 0.0f);
    if (std::abs(pxr::GfDot(tangent, axis)) > 0.98f)
    {
        axis = pxr::GfVec3f(1.0f, 0.0f, 0.0f);
    }
    pxr::GfVec3f n = pxr::GfCross(tangent, axis);
    float len = n.GetLength();
    if (len <= 1e-8f)
    {
        axis = pxr::GfVec3f(0.0f, 0.0f, 1.0f);
        n = pxr::GfCross(tangent, axis);
        len = n.GetLength();
    }
    return (len > 1e-8f) ? (n / len) : pxr::GfVec3f(1.0f, 0.0f, 0.0f);
}
} // namespace

bool WriteEdgeRateDebugObj(const std::string &path,
                           const std::vector<EdgeRateDebugLine> &lines,
                           int minRate,
                           int maxRate)
{
    std::ofstream out(path);
    if (!out.is_open())
    {
        return false;
    }

    out << "# edge_rate_debug\n";
    out << "# triangles generated from line segments for viewport color preview\n";
    const float halfWidth = 0.5f * ComputeDebugWidth(lines);
    int nextIndex = 1;
    for (const EdgeRateDebugLine &line : lines)
    {
        if (line.points.size() < 2)
        {
            continue;
        }
        const pxr::GfVec3f c = RateToColor(line.rate, minRate, maxRate);
        for (int i = 0; i + 1 < int(line.points.size()); ++i)
        {
            const pxr::GfVec3f p0 = line.points[size_t(i)];
            const pxr::GfVec3f p1 = line.points[size_t(i + 1)];
            const pxr::GfVec3f d = p1 - p0;
            const float len = d.GetLength();
            if (len <= 1e-8f)
            {
                continue;
            }
            const pxr::GfVec3f t = d / len;
            const pxr::GfVec3f n = SafePerp(t) * halfWidth;

            const pxr::GfVec3f v0 = p0 - n;
            const pxr::GfVec3f v1 = p0 + n;
            const pxr::GfVec3f v2 = p1 + n;
            const pxr::GfVec3f v3 = p1 - n;

            out << "v " << v0[0] << " " << v0[1] << " " << v0[2] << " " << c[0] << " " << c[1] << " "
                << c[2] << "\n";
            out << "v " << v1[0] << " " << v1[1] << " " << v1[2] << " " << c[0] << " " << c[1] << " "
                << c[2] << "\n";
            out << "v " << v2[0] << " " << v2[1] << " " << v2[2] << " " << c[0] << " " << c[1] << " "
                << c[2] << "\n";
            out << "v " << v3[0] << " " << v3[1] << " " << v3[2] << " " << c[0] << " " << c[1] << " "
                << c[2] << "\n";
            out << "f " << nextIndex << " " << (nextIndex + 1) << " " << (nextIndex + 2) << "\n";
            out << "f " << nextIndex << " " << (nextIndex + 2) << " " << (nextIndex + 3) << "\n";
            nextIndex += 4;
        }
    }
    return out.good();
}
} // namespace ybi::testio
