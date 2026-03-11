#include "util/math_common.h"
#include "util/vec3.h"

YBI_NAMESPACE_BEGIN

static bool ComputeTriangleShadingNormal(const LaunchParams &params,
                                         float u,
                                         float v,
                                         int primitiveIndex,
                                         int instanceId,
                                         Vec3 &outNormal)
{
    assert(outNormal);
    const LaunchParams::InstanceGeomRef *refs =
        reinterpret_cast<const LaunchParams::InstanceGeomRef *>(params.instanceGeomRefs);
    const LaunchParams::InstanceGeomRef ref = refs[instanceId];

    const int triCornerBase = primitiveIndex * 3;
    assert(triCornerBase + 2 < ref.numNormalIndices);

    const Vec3 *normals = reinterpret_cast<const Vec3 *>(ref.normals);
    const int *normalIndices = reinterpret_cast<const int *>(ref.normalIndices);
    const int n0 = normalIndices[triCornerBase + 0];
    const int n1 = normalIndices[triCornerBase + 1];
    const int n2 = normalIndices[triCornerBase + 2];
    assert(!(n0 < 0 || n0 >= ref.numNormals || n1 < 0 || n1 >= ref.numNormals || n2 < 0 ||
             n2 >= ref.numNormals));

    const float u = rayHit.hit.u;
    const float v = rayHit.hit.v;
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

YBI_NAMESPACE_END
