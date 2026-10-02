#pragma once

#include "core/affine_transform.h"
#include "core/frame.h"
#include "core/vec2f.h"
#include "core/vec3f.h"
#include "render/bsdf/materialx/bsdf_types.h"
#include "render/types.h"

namespace ybi
{

struct ShadingInfo
{
    Vec3f positions[3];
    Vec3f normals[3];
    Vec2f uvs[3];

    Vec3f position;
    Vec3f geometricNormal;
    Vec3f shadingNormal;
    Vec3f dpdu, dpdv;
    Vec3f dndu, dndv;
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

    float closureData[MAX_CLOSURE_SIZE];

    YBI_GPU void TransformToWorld(const AffineTransform &objectToWorld)
    {
        position = objectToWorld.TransformPoint(position);
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
            if (!hasNormals)
                continue;
            normals[i] = Normalize(Vec3f(Dot(adjugate0, normals[i]),
                                         Dot(adjugate1, normals[i]),
                                         Dot(adjugate2, normals[i])));
        }

        shadingNormal = Normalize(Vec3f(Dot(adjugate0, shadingNormal),
                                        Dot(adjugate1, shadingNormal),
                                        Dot(adjugate2, shadingNormal)));
        geometricNormal = Normalize(Vec3f(Dot(adjugate0, geometricNormal),
                                          Dot(adjugate1, geometricNormal),
                                          Dot(adjugate2, geometricNormal)));

        dpdu = objectToWorld.TransformVector(dpdu);
        dpdv = objectToWorld.TransformVector(dpdv);

        const float invDet = 1.f / Dot(r0, adjugate0);
        dndu = Vec3f(Dot(adjugate0, dndu), Dot(adjugate1, dndu), Dot(adjugate2, dndu)) * invDet;
        dndv = Vec3f(Dot(adjugate0, dndv), Dot(adjugate1, dndv), Dot(adjugate2, dndv)) * invDet;
    }
};

template <typename T>
YBI_GPU static T InterpolateTriangleValues(T values[3], Vec2f bary)
{
    return values[0] * (1.f - bary.x - bary.y) + values[1] * bary.x + values[2] * bary.y;
}

YBI_CPU_GPU void CalculateTriangleHitInfo(ShadingInfo &shadingInfo)
{
    const Vec3f &p0 = shadingInfo.positions[0];
    const Vec3f &p1 = shadingInfo.positions[1];
    const Vec3f &p2 = shadingInfo.positions[2];
    const Vec3f gn = Normalize(Cross(p1 - p0, p2 - p0));
    const Vec3f n0 = shadingInfo.hasNormals ? shadingInfo.normals[0] : gn;
    const Vec3f n1 = shadingInfo.hasNormals ? shadingInfo.normals[1] : gn;
    const Vec3f n2 = shadingInfo.hasNormals ? shadingInfo.normals[2] : gn;
    const Vec2f uv0 = shadingInfo.hasUvs ? shadingInfo.uvs[0] : Vec2f(0.f);
    const Vec2f uv1 = shadingInfo.hasUvs ? shadingInfo.uvs[1] : Vec2f(0.f);
    const Vec2f uv2 = shadingInfo.hasUvs ? shadingInfo.uvs[2] : Vec2f(0.f);
    const Vec2f &bary = shadingInfo.bary;

    Vec2f duv10 = uv1 - uv0;
    Vec2f duv20 = uv2 - uv0;
    float det = fmaf(duv10.x, duv20.y, -duv10.y * duv20.x);

    Vec3f dpdu(0.f), dpdv(0.f), dndu(0.f), dndv(0.f);
    Vec3f dp10 = p1 - p0;
    Vec3f dp20 = p2 - p0;
    Vec3f dn10 = n1 - n0;
    Vec3f dn20 = n2 - n0;

    Vec3f origin = p0 + dp10 * bary.x + dp20 * bary.y;
    Vec2f uv = uv0 + duv10 * bary.x + duv20 * bary.y;

    if (Abs(det) < 1e-9f)
    {
        BuildOrthonormalBasis(gn, dpdu, dpdv);
    }
    else
    {
        float invDet = 1.f / det;

        dpdu = (duv20.y * dp10 - duv10.y * dp20) * invDet;
        dpdv = (-duv20.x * dp10 + duv10.x * dp20) * invDet;

        dndu = (duv20.y * dn10 - duv10.y * dn20) * invDet;
        dndv = (-duv20.x * dn10 + duv10.x * dn20) * invDet;
    }

    Vec3f n = Normalize(n0 + dn10 * bary[0] + dn20 * bary[1]);

    Vec3f ss = dpdu;
    Vec3f ts = Cross(n, ss);
    if (Dot(ts, ts) > 0)
    {
        ss = Cross(ts, n);
    }
    else
    {
        BuildOrthonormalBasis(n, ss, ts);
    }

    shadingInfo.position = origin;
    shadingInfo.uv = uv;
    shadingInfo.shadingNormal = n;
    shadingInfo.geometricNormal = gn;
    shadingInfo.dpdu = ss;
    shadingInfo.dpdv = ts;
    shadingInfo.dndu = dndu;
    shadingInfo.dndv = dndv;
}

} // namespace ybi
