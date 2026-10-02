#pragma once

#include "kernel/geometry/clustered_triangle.h"
#include "kernel/geometry/triangle.h"
#include "kernel/integrator/ray_cone.h"
#include "kernel/shading/shading.h"
#include "kernel/launch/shading_context.h"
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
void ShadeHit(const ShadingContext &context, SampleState *sampleState, uint32_t sampleIndex)
{
    const DeviceScene &scene = context.scene;
    Vec3f rayDirection = sampleState->ReadRayDirection(sampleIndex);
    Vec3f rayOrigin = sampleState->ReadRayOrigin(sampleIndex);
    RayCone rayCone = sampleState->ReadRayCone(sampleIndex);
    const IntersectionInfo info = sampleState->ReadIntersection(sampleIndex);
    uint32_t depth = sampleState->ReadDepth(sampleIndex);

    const Object &object = scene.objects[info.id];

    ShadingInfo shadingInfo;
    GetShadingInfo(scene, object, info, shadingInfo);
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

    const Vec3f worldViewDirection = -rayDirection;
    const DeviceMaterial3 *materialXMaterial = nullptr;

    bool passedAlphaTest = true;
    bool hasBsdf = false;

    if (shadingInfo.hasMaterial)
    {
        const MaterialHandle &handle = scene.materialHandles[shadingInfo.materialIndex];
        assert(handle.type == MaterialHandleType::MaterialX);
        materialXMaterial = &scene.deviceMaterial3s[handle.index];

        DeviceDispatch(materialXMaterial->initializeCallableIndex, &params, &shadingInfo);
        hasBsdf = true;
    }

    if (!shadingInfo.hasMaterial || !hasBsdf || !passedAlphaTest)
    {
        rayOrigin = OffsetRayOrigin(shadingInfo.position, -geometricNormal);
        depth++;
        // return;
    }
    else
    {
    }

    sampleState->WriteRayOrigin(sampleIndex, rayOrigin);
    sampleState->WriteRayCone(sampleIndex, rayCone);
    sampleState->WriteDepth(sampleIndex, depth);
}

} // namespace ybi
