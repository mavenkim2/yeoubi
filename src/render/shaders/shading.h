#pragma once

#include "core/affine_transform.h"
#include "core/vec2f.h"
#include "core/vec3f.h"

namespace ybi
{

struct ShadingInfo
{
    Vec3f positions[3];
    Vec3f normals[3];
    Vec2f uvs[3];

    Vec3f shadingNormal;
    Vec2f uv;
    Vec2f bary;

    Vec3f rayDirection;
    RayCone rayCone;

    uint32_t objectId;
    uint32_t clusterId;
    uint32_t primIndex;

    uint32_t materialIndex;
    bool hasNormals;
    bool hasUvs;
    bool hasMaterial;

    YBI_GPU void TransformToWorld(const AffineTransform &objectToWorld)
    {
        for (int i = 0; i < 3; i++)
        {
            positions[i] = objectToWorld.TransformPoint(positions[i]);
        }

        Vec3f r0 = objectToWorld.GetRow(0).xyz;
        Vec3f r1 = objectToWorld.GetRow(1).xyz;
        Vec3f r2 = objectToWorld.GetRow(2).xyz;

        Vec3f adjugate0 = Cross(r1, r2);
        Vec3f adjugate1 = Cross(r2, r0);
        Vec3f adjugate2 = Cross(r0, r1);

        for (int i = 0; i < 3; i++)
        {
            normals[i] = Normalize(Vec3f(Dot(adjugate0, normals[i]),
                                         Dot(adjugate1, normals[i]),
                                         Dot(adjugate2, normals[i])));
        }

        shadingNormal = Normalize(Vec3f(Dot(adjugate0, shadingNormal),
                                        Dot(adjugate1, shadingNormal),
                                        Dot(adjugate2, shadingNormal)));
    }
};

template <typename T>
YBI_GPU static T InterpolateTriangleValues(T values[3], Vec2f bary)
{
    return values[0] * (1.f - bary.x - bary.y) + values[1] * bary.x + values[2] * bary.y;
}

} // namespace ybi
