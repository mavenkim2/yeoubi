#pragma once

#include "kernel/backend/dispatch.h"
#include "kernel/geometry/clustered_triangle.h"
#include "kernel/geometry/triangle.h"
#include "kernel/launch/kernel_params.h"
#include "kernel/raytracing/ray_offset.h"
#include "kernel/shading/shading.h"
#include "kernel/state/types.h"

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

template <typename SampleState>
YBI_GPU void
ShadeHit(const KernelParams &params, const SampleState *sampleState, uint32_t sampleIndex)
{
    const DeviceScene &scene = params.scene;
    const Vec3f rayDirection = sampleState->ReadRayDirection(sampleIndex);
    Vec3f rayOrigin = sampleState->ReadRayOrigin(sampleIndex);
    uint32_t depth = sampleState->ReadDepth(sampleIndex);
    const Vec3f position = sampleState->ReadPosition(sampleIndex);
    const uint32_t materialIndex = sampleState->ReadMaterialIndex(sampleIndex);
    Vec3f geometricNormal = sampleState->ReadGeometricNormal(sampleIndex);
    ClosureInitializationInfo closureInitializationInfo =
        sampleState->ReadClosureInitializationInfo(sampleIndex);

    if (Dot(geometricNormal, rayDirection) > 0.f)
        geometricNormal = -geometricNormal;

    bool passedAlphaTest = true;
    bool hasBsdf = false;
    const bool hasMaterial = materialIndex != ~0u;

    if (hasMaterial)
    {
        const MaterialHandle &handle = scene.materialHandles[materialIndex];
        assert(handle.type == MaterialHandleType::MaterialX);
        const DeviceMaterial3 &materialXMaterial = scene.deviceMaterial3s[handle.index];

        DeviceDispatch<void>(
            materialXMaterial.initializeCallableIndex, &params, &closureInitializationInfo);
        hasBsdf = true;
    }

    if (!hasMaterial || !hasBsdf || !passedAlphaTest)
    {
        rayOrigin = OffsetRayOrigin(position, -geometricNormal);
        depth++;
        sampleState->WriteKernelType(sampleIndex, KERNEL_INTERSECT);
        sampleState->WriteRayOrigin(sampleIndex, rayOrigin);
        sampleState->WriteDepth(sampleIndex, depth);
        // return;
    }
    else
    {
        sampleState->WriteKernelType(sampleIndex, KERNEL_DIRECT_LIGHTING);
    }
}

} // namespace ybi
