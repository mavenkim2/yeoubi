#pragma once

#include "core/affine_transform.h"
#include "core/vec2f.h"
#include "render/ray_cone.h"
#include "render/rng.h"

namespace ybi
{

enum IntersectionType
{
    INTERSECTION_TYPE_NONE = 0,
    INTERSECTION_TYPE_TRIANGLE = 1,
};

struct IntersectionInfo
{
    float t;
    Vec2f bary;
    AffineTransform objectToWorld;
    uint32_t primIndex;
    uint32_t clusterId;
    uint32_t id;
    IntersectionType type;
};

struct IntersectionSampleState
{
    IntersectionInfo info;
    uint32_t sampleIndex;
};

struct SampleState
{
    uint32_t pixelX;
    uint32_t pixelY;
    RNG rng;
    RayCone rayCone;

    Vec3f radiance;
    Vec3f throughput;
    Vec3f previousSurfacePosition;
    float previousBsdfPdf;
};

} // namespace ybi
