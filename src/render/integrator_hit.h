#pragma once

#include "render/integrator_common.h"
#include "render/launch_params.h"
#include "util/math_common.h"
#include "util/vec3.h"
#include <assert.h>

namespace ybi
{

static YBI_INTEGRATOR_HD bool
TryGetTriangleLocalPositions(const LaunchParams &params,
                             int instanceId,
                             int primitiveIndex,
                             Vec3 &outP0,
                             Vec3 &outP1,
                             Vec3 &outP2)
{
    if (params.instanceGeomRefs == 0ull || instanceId < 0 ||
        instanceId >= params.instanceGeomRefCount || primitiveIndex < 0)
    {
        return false;
    }

    const LaunchParams::InstanceGeomRef *refs =
        reinterpret_cast<const LaunchParams::InstanceGeomRef *>(params.instanceGeomRefs);
    const LaunchParams::InstanceGeomRef ref = refs[instanceId];
    if (ref.positions == 0ull || ref.indices == 0ull)
    {
        return false;
    }

    const int triCornerBase = primitiveIndex * 3;
    if (triCornerBase + 2 >= ref.numIndices)
    {
        return false;
    }

    const Vec3 *positions = reinterpret_cast<const Vec3 *>(ref.positions);
    const int *indices = reinterpret_cast<const int *>(ref.indices);
    const int i0 = indices[triCornerBase + 0];
    const int i1 = indices[triCornerBase + 1];
    const int i2 = indices[triCornerBase + 2];
    if (i0 < 0 || i0 >= ref.numPositions || i1 < 0 || i1 >= ref.numPositions || i2 < 0 ||
        i2 >= ref.numPositions)
    {
        return false;
    }

    outP0 = positions[i0];
    outP1 = positions[i1];
    outP2 = positions[i2];
    return true;
}

static YBI_INTEGRATOR_HD bool ComputeTriangleShadingNormal(const LaunchParams &params,
                                                           float u,
                                                           float v,
                                                           int primitiveIndex,
                                                           int instanceId,
                                                           Vec3 &outNormal)
{
    const LaunchParams::InstanceGeomRef *refs =
        reinterpret_cast<const LaunchParams::InstanceGeomRef *>(params.instanceGeomRefs);
    const LaunchParams::InstanceGeomRef ref = refs[instanceId];

    if (ref.numNormalIndices == 0)
        return false;

    const int triCornerBase = primitiveIndex * 3;
    assert(triCornerBase + 2 < ref.numNormalIndices);

    const Vec3 *normals = reinterpret_cast<const Vec3 *>(ref.normals);
    const int *normalIndices = reinterpret_cast<const int *>(ref.normalIndices);
    const int n0 = normalIndices[triCornerBase + 0];
    const int n1 = normalIndices[triCornerBase + 1];
    const int n2 = normalIndices[triCornerBase + 2];
    assert(!(n0 < 0 || n0 >= ref.numNormals || n1 < 0 || n1 >= ref.numNormals || n2 < 0 ||
             n2 >= ref.numNormals));

    const float w = 1.0f - u - v;
    const Vec3 localNormal(normals[n0].x * w + normals[n1].x * u + normals[n2].x * v,
                           normals[n0].y * w + normals[n1].y * u + normals[n2].y * v,
                           normals[n0].z * w + normals[n1].z * u + normals[n2].z * v);
    if (Length(localNormal) <= 1.0e-12f)
    {
        return false;
    }

    outNormal = localNormal;
    return true;
}

static YBI_INTEGRATOR_HD bool TryGetTriangleTexcoords(
    const LaunchParams &params,
    int instanceId,
    int primitiveIndex,
    render::integrator::UV2 &outUv0,
    render::integrator::UV2 &outUv1,
    render::integrator::UV2 &outUv2)
{
    if (params.instanceGeomRefs == 0ull || instanceId < 0 ||
        instanceId >= params.instanceGeomRefCount || primitiveIndex < 0)
    {
        return false;
    }

    const LaunchParams::InstanceGeomRef *refs =
        reinterpret_cast<const LaunchParams::InstanceGeomRef *>(params.instanceGeomRefs);
    const LaunchParams::InstanceGeomRef ref = refs[instanceId];
    if (ref.texcoords == 0ull || ref.texcoordIndices == 0ull)
    {
        return false;
    }

    const int triCornerBase = primitiveIndex * 3;
    if (triCornerBase + 2 >= ref.numTexcoordIndices)
    {
        return false;
    }

    const int *texcoordIndices = reinterpret_cast<const int *>(ref.texcoordIndices);
    const int t0 = texcoordIndices[triCornerBase + 0];
    const int t1 = texcoordIndices[triCornerBase + 1];
    const int t2 = texcoordIndices[triCornerBase + 2];
    if (t0 < 0 || t0 >= ref.numTexcoords || t1 < 0 || t1 >= ref.numTexcoords || t2 < 0 ||
        t2 >= ref.numTexcoords)
    {
        return false;
    }

    const render::integrator::UV2 *texcoords =
        reinterpret_cast<const render::integrator::UV2 *>(ref.texcoords);
    outUv0 = texcoords[t0];
    outUv1 = texcoords[t1];
    outUv2 = texcoords[t2];
    return true;
}

static YBI_INTEGRATOR_HD bool TryComputeTriangleSurfacePartials(
    const LaunchParams &params, render::integrator::HitInfo *outHit)
{
    if (!outHit || !outHit->hasWorldTriangle || outHit->instanceId < 0 || outHit->primitiveIndex < 0)
    {
        return false;
    }

    outHit->dPds = Vec3(0.0f, 0.0f, 0.0f);
    outHit->dPdt = Vec3(0.0f, 0.0f, 0.0f);
    outHit->hasSurfacePartials = false;

    render::integrator::UV2 uv0 = {};
    render::integrator::UV2 uv1 = {};
    render::integrator::UV2 uv2 = {};
    if (!TryGetTriangleTexcoords(
            params, outHit->instanceId, outHit->primitiveIndex, uv0, uv1, uv2))
    {
        return false;
    }

    const Vec3 edge1 = outHit->worldTri1 - outHit->worldTri0;
    const Vec3 edge2 = outHit->worldTri2 - outHit->worldTri0;
    const float du1 = uv1.x - uv0.x;
    const float dv1 = uv1.y - uv0.y;
    const float du2 = uv2.x - uv0.x;
    const float dv2 = uv2.y - uv0.y;
    const float det = du1 * dv2 - dv1 * du2;
    if (fabsf(det) <= 1.0e-8f)
    {
        return false;
    }

    const float invDet = 1.0f / det;
    outHit->dPds = (edge1 * dv2 - edge2 * dv1) * invDet;
    outHit->dPdt = (edge2 * du1 - edge1 * du2) * invDet;
    outHit->hasSurfacePartials = true;
    return true;
}

} // namespace ybi
