#pragma once

#include "render/kernels/clustered_triangle.h"
#include "render/kernels/triangle.h"
#include "render/shaders/shading.h"

namespace ybi
{

YBI_GPU void GetShadingInfo(const DeviceScene &deviceScene,
                            const Object &object,
                            const IntersectionInfo &info,
                            ShadingInfo &shadingInfo)
{
    shadingInfo = {};
    if (IsClusteredTriangle(object.clusterOffset, info.clusterId))
    {
        const uint32_t clusterIndex = object.clusterOffset + info.clusterId;
        const DeviceCluster &cluster = deviceScene.compressedClusters[clusterIndex];
        const uint8_t *clusterData = deviceScene.compressedClusterBytes + cluster.byteOffset;
        const Vec3u vertexIndices = cluster.DecodeTriangle(clusterData, info.primIndex);
        GetClusteredShadingInfo(object, info, cluster, clusterData, vertexIndices, shadingInfo);
    }
    else
    {
        GetTriangleShadingInfo(deviceScene, object, info, shadingInfo);
    }

    shadingInfo.bary = info.bary;
    CalculateTriangleHitInfo(shadingInfo);
    shadingInfo.objectId = info.id;
    shadingInfo.clusterId = info.clusterId;
    shadingInfo.primIndex = info.primIndex;
}

#if 0
void ShadeHit(const DeviceScene &scene)
{
    const Object &object = params.scene.objects[info.id];

    ShadingInfo shadingInfo;
    GetShadingInfo(object, info, shadingInfo);
    shadingInfo.TransformToWorld(info.objectToWorld);

    Vec3f geometricNormal = shadingInfo.geometricNormal;
    if (Dot(geometricNormal, rayDirection) > 0.f)
        geometricNormal = -geometricNormal;

    // if (Dot(geometricNormal, rayDirection) > 0.f)
    // {
    //     shadingInfo.shadingNormal = -shadingInfo.shadingNormal;
    // }

    UpdateRayCone(rayCone, rayDirection, info.t, shadingInfo, depth);
    shadingInfo.rayDirection = rayDirection;
    shadingInfo.rayCone = rayCone;

#if defined(YBI_MATERIAL_BACKEND_OPENPBR)
    OpenPbrBsdf bsdf;
#endif

#if defined(YBI_MATERIAL_BACKEND_OPENPBR)
    const OrthonormalBasis tangentFrame(shadingInfo.dpdu, shadingInfo.shadingNormal);
#endif
    const Vec3f worldViewDirection = -rayDirection;
#if defined(YBI_MATERIAL_BACKEND_MATERIALX)
    const DeviceMaterial3 *materialXMaterial = nullptr;
#endif

    bool passedAlphaTest = true;
    bool hasBsdf = false;

    if (shadingInfo.hasMaterial)
    {
        const MaterialHandle &handle = params.scene.materialHandles[shadingInfo.materialIndex];
#if defined(YBI_MATERIAL_BACKEND_OPENPBR)
        if (handle.type == MaterialHandleType::OpenPbr)
        {
            const DeviceMaterial2 &material = params.scene.deviceMaterial2s[handle.index];

            passedAlphaTest = InitializeOpenPbrBsdf(params,
                                                    material,
                                                    shadingInfo,
                                                    bsdf,
                                                    shadingInfo.shadingNormal,
                                                    worldViewDirection,
                                                    tangentFrame,
                                                    throughput,
                                                    rng.Uniform());
            hasBsdf = true;
        }
#else
        if (handle.type == MaterialHandleType::MaterialX)
        {
            materialXMaterial = &params.scene.deviceMaterial3s[handle.index];
            optixDirectCall<void>(
                materialXMaterial->initializeCallableIndex, &params, &shadingInfo);
            hasBsdf = true;
        }
#endif
    }

    if (!shadingInfo.hasMaterial || !hasBsdf || !passedAlphaTest)
    {
        rayOrigin = OffsetRayOrigin(shadingInfo.position, -geometricNormal);
        ++depth;
        continue;
    }
}
#endif

} // namespace ybi
