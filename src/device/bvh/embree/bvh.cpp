#include "device/cpu_device.h"

#if defined(WITH_EMBREE)

#include "scene/scene.h"
#include "util/assert.h"

#include <cstdio>
#include <cstring>

YBI_NAMESPACE_BEGIN

namespace
{
static void CopyTransformMatrix(const float3x4 &src, float dst[12])
{
    std::memcpy(dst, src.m, sizeof(float) * 12);
}

static RTCScene BuildTriangleLeafScene(RTCDevice embreeDevice, const Mesh &mesh)
{
    YBI_ASSERT(embreeDevice);
    YBI_ASSERT(mesh.positions.size() > 0);
    YBI_ASSERT(mesh.indices.size() >= 3);

    RTCScene leafScene = rtcNewScene(embreeDevice);
    YBI_ASSERT(leafScene);

    RTCGeometry triangles = rtcNewGeometry(embreeDevice, RTC_GEOMETRY_TYPE_TRIANGLE);
    YBI_ASSERT(triangles);

    rtcSetSharedGeometryBuffer(triangles,
                               RTC_BUFFER_TYPE_VERTEX,
                               0,
                               RTC_FORMAT_FLOAT3,
                               mesh.positions.data(),
                               0,
                               sizeof(float3),
                               mesh.positions.size());
    rtcSetSharedGeometryBuffer(triangles,
                               RTC_BUFFER_TYPE_INDEX,
                               0,
                               RTC_FORMAT_UINT3,
                               mesh.indices.data(),
                               0,
                               sizeof(int) * 3,
                               mesh.indices.size() / 3);

    rtcCommitGeometry(triangles);
    rtcAttachGeometry(leafScene, triangles);
    rtcReleaseGeometry(triangles);
    rtcCommitScene(leafScene);
    return leafScene;
}
} // namespace

void BuildEmbreeBVH(CPUDevice *cpuDevice, Scene *scene)
{
    YBI_ASSERT(cpuDevice);
    YBI_ASSERT(scene);
    YBI_ASSERT(cpuDevice->embreeDevice);

    if (scene->bvhHandle)
    {
        rtcReleaseScene((RTCScene)scene->bvhHandle);
        scene->bvhHandle = 0;
    }

    RTCScene embreeScene = rtcNewScene(cpuDevice->embreeDevice);
    YBI_ASSERT(embreeScene);

    for (const Mesh &mesh : scene->meshes)
    {
        if (mesh.positions.size() == 0 || mesh.indices.size() == 0)
        {
            continue;
        }

        YBI_ASSERT(mesh.refIndex != UINT32_MAX);

        RTCScene leafScene = BuildTriangleLeafScene(cpuDevice->embreeDevice, mesh);

        RTCGeometry instance = rtcNewGeometry(cpuDevice->embreeDevice, RTC_GEOMETRY_TYPE_INSTANCE);
        YBI_ASSERT(instance);
        rtcSetGeometryInstancedScene(instance, leafScene);
        float transform[12] = {};
        CopyTransformMatrix(mesh.parentFromLocal, transform);
        rtcSetGeometryTransform(instance, 0, RTC_FORMAT_FLOAT3X4_ROW_MAJOR, transform);
        rtcCommitGeometry(instance);
        rtcAttachGeometryByID(embreeScene, instance, mesh.refIndex);
        rtcReleaseGeometry(instance);

        rtcReleaseScene(leafScene);
    }

    static bool warnedCurves = false;
    for (const Curves &curves : scene->curves)
    {
        if (curves.GetNumVertices() == 0 || curves.GetNumCurves() == 0)
        {
            continue;
        }
        if (!warnedCurves)
        {
            fprintf(stderr, "Embree CPU path: curve geometry build TODO; skipping curves.\n");
            warnedCurves = true;
        }
    }

    static bool warnedMissingChild = false;
    for (const Instance &instanceData : scene->instances)
    {
        YBI_ASSERT(instanceData.childSceneIndex < scene->childScenes.size());
        Scene *childScene = scene->childScenes[instanceData.childSceneIndex];
        YBI_ASSERT(childScene);
        RTCScene childEmbreeScene = (RTCScene)childScene->bvhHandle;
        if (!childEmbreeScene)
        {
            if (!warnedMissingChild)
            {
                fprintf(stderr,
                        "Embree CPU path: child scene BVH missing; skipping instance.\n");
                warnedMissingChild = true;
            }
            continue;
        }

        RTCGeometry childInstance =
            rtcNewGeometry(cpuDevice->embreeDevice, RTC_GEOMETRY_TYPE_INSTANCE);
        YBI_ASSERT(childInstance);
        rtcSetGeometryInstancedScene(childInstance, childEmbreeScene);
        float transform[12] = {};
        CopyTransformMatrix(instanceData.parentFromLocal, transform);
        rtcSetGeometryTransform(childInstance, 0, RTC_FORMAT_FLOAT3X4_ROW_MAJOR, transform);
        rtcCommitGeometry(childInstance);
        rtcAttachGeometry(embreeScene, childInstance);
        rtcReleaseGeometry(childInstance);
    }

    rtcCommitScene(embreeScene);
    scene->bvhHandle = (Scene::BVHHandle)embreeScene;
}

YBI_NAMESPACE_END

#endif

