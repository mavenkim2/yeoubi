#pragma once

#include "core/affine_transform.h"
#include "core/frame.h"
#include "core/math.h"
#include "kernel/geometry/clustered_triangle.h"
#include "kernel/geometry/triangle.h"
#include "kernel/state/types.h"

namespace ybi
{

struct TriangleInfo
{
    Vec3f positions[3];
    Vec3f normals[3];
    Vec2f uvs[3];

    uint32_t materialIndex;
    bool hasNormals;
    bool hasUvs;
    bool hasMaterial;

    YBI_GPU void TransformToWorld(const AffineTransform &objectToWorld)
    {
        for (uint32_t i = 0; i < 3u; ++i)
        {
            positions[i] = objectToWorld.TransformPoint(positions[i]);
        }

        if (!hasNormals)
        {
            return;
        }

        const Vec3f r0 = objectToWorld.GetRow(0).xyz;
        const Vec3f r1 = objectToWorld.GetRow(1).xyz;
        const Vec3f r2 = objectToWorld.GetRow(2).xyz;
        const Vec3f adjugate0 = Cross(r1, r2);
        const Vec3f adjugate1 = Cross(r2, r0);
        const Vec3f adjugate2 = Cross(r0, r1);

        for (uint32_t i = 0; i < 3u; ++i)
        {
            normals[i] = Normalize(Vec3f(Dot(adjugate0, normals[i]),
                                         Dot(adjugate1, normals[i]),
                                         Dot(adjugate2, normals[i])));
        }
    }
};

YBI_GPU void GetTriangleInfo(const DeviceScene &scene,
                             const Object &object,
                             const IntersectionInfo &intersection,
                             TriangleInfo &triangleInfo)
{
    triangleInfo = {};
    triangleInfo.materialIndex = ~0u;

    if (IsClusteredTriangle(object.clusterOffset, intersection.clusterId))
    {
        const uint32_t clusterIndex = object.clusterOffset + intersection.clusterId;
        const DeviceCluster &cluster = scene.compressedClusters[clusterIndex];
        const uint8_t *clusterData = scene.compressedClusterBytes + cluster.byteOffset;
        const Vec3u vertexIndices =
            cluster.DecodeTriangle(clusterData, intersection.primIndex);

        GetClusteredTrianglePositions(
            cluster, clusterData, vertexIndices, triangleInfo.positions);
        triangleInfo.hasNormals = GetClusteredTriangleNormals(
            object, intersection, cluster, clusterData, vertexIndices, triangleInfo.normals);
        triangleInfo.hasUvs = GetClusteredTriangleUVs(
            object, intersection, cluster, clusterData, vertexIndices, triangleInfo.uvs);

        if (object.materialIndex != -1)
        {
            triangleInfo.materialIndex = object.materialIndex;
            triangleInfo.hasMaterial = true;
        }
        else
        {
            triangleInfo.hasMaterial = GetClusteredTriangleMaterial(
                intersection, cluster, clusterData, triangleInfo.materialIndex);
        }
        return;
    }

    const Vec3u vertexIndices = GetTriangleVertexIndices(scene, object, intersection);
    GetTrianglePositions(scene, object, vertexIndices, triangleInfo.positions);
    triangleInfo.hasNormals =
        GetTriangleNormals(scene, object, intersection, vertexIndices, triangleInfo.normals);
    triangleInfo.hasUvs =
        GetTriangleUVs(scene, object, intersection, vertexIndices, triangleInfo.uvs);
    triangleInfo.hasMaterial =
        GetTriangleMaterial(scene, object, intersection.primIndex, triangleInfo.materialIndex);
}

YBI_GPU void CalculateTriangleSurfaceInfo(const TriangleInfo &triangleInfo,
                                          const Vec2f &bary,
                                          TriangleSurfaceInfo &surfaceInfo)
{
    const Vec3f &p0 = triangleInfo.positions[0];
    const Vec3f &p1 = triangleInfo.positions[1];
    const Vec3f &p2 = triangleInfo.positions[2];
    const Vec3f dp10 = p1 - p0;
    const Vec3f dp20 = p2 - p0;
    const Vec3f geometricNormalUnnormalized = Cross(dp10, dp20);
    const Vec3f geometricNormal = Normalize(geometricNormalUnnormalized);

    const Vec3f n0 = triangleInfo.hasNormals ? triangleInfo.normals[0] : geometricNormal;
    const Vec3f n1 = triangleInfo.hasNormals ? triangleInfo.normals[1] : geometricNormal;
    const Vec3f n2 = triangleInfo.hasNormals ? triangleInfo.normals[2] : geometricNormal;
    const Vec3f dn10 = n1 - n0;
    const Vec3f dn20 = n2 - n0;

    const Vec2f uv0 = triangleInfo.hasUvs ? triangleInfo.uvs[0] : Vec2f(0.f);
    const Vec2f uv1 = triangleInfo.hasUvs ? triangleInfo.uvs[1] : Vec2f(0.f);
    const Vec2f uv2 = triangleInfo.hasUvs ? triangleInfo.uvs[2] : Vec2f(0.f);
    const Vec2f duv10 = uv1 - uv0;
    const Vec2f duv20 = uv2 - uv0;
    const float uvDeterminant = fmaf(duv10.x, duv20.y, -duv10.y * duv20.x);

    Vec3f dpdu;
    if (Abs(uvDeterminant) < 1e-9f)
    {
        Vec3f dpdv;
        BuildOrthonormalBasis(geometricNormal, dpdu, dpdv);
    }
    else
    {
        const float inverseUvDeterminant = 1.f / uvDeterminant;
        dpdu = (duv20.y * dp10 - duv10.y * dp20) * inverseUvDeterminant;
    }

    const Vec3f shadingNormal = Normalize(n0 + dn10 * bary.x + dn20 * bary.y);
    Vec3f tangent = dpdu;
    Vec3f bitangent = Cross(shadingNormal, tangent);
    if (Dot(bitangent, bitangent) > 0.f)
    {
        tangent = Cross(bitangent, shadingNormal);
    }
    else
    {
        BuildOrthonormalBasis(shadingNormal, tangent, bitangent);
    }

    surfaceInfo.position = p0 + dp10 * bary.x + dp20 * bary.y;
    surfaceInfo.geometricNormal = geometricNormal;
    surfaceInfo.shadingNormal = shadingNormal;
    surfaceInfo.dpdu = tangent;
    surfaceInfo.uv = uv0 + duv10 * bary.x + duv20 * bary.y;
    surfaceInfo.triangleArea = Length(geometricNormalUnnormalized);
    surfaceInfo.uvArea = Abs(uvDeterminant);
    surfaceInfo.materialIndex = triangleInfo.hasMaterial ? triangleInfo.materialIndex : ~0u;
}

} // namespace ybi
