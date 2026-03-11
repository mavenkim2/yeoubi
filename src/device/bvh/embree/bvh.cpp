#include "device/cpu_device.h"

#if defined(WITH_EMBREE)

#include "scene/scene.h"
#include "util/assert.h"
#include "util/vec4.h"

#include <algorithm>
#include <cstdint>
#include <cstdio>
#include <cstring>

namespace ybi
{

namespace
{
static void CopyTransformMatrix(const Float3x4 &src, float dst[12])
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
                               sizeof(Vec3),
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

static float ResolveCurveWidth(const Curves &curves,
                               const Array<float> &widths,
                               uint32_t curveIndex,
                               uint32_t vertexIndex)
{
    if (widths.size() == 0)
    {
        return 0.01f;
    }
    if (widths.size() == curves.GetNumVertices())
    {
        return widths[vertexIndex];
    }
    if (widths.size() == curves.GetNumCurves())
    {
        return widths[curveIndex];
    }
    return widths[0];
}

static RTCScene BuildCurveLeafScene(RTCDevice embreeDevice, const Curves &curves)
{
    YBI_ASSERT(embreeDevice);
    YBI_ASSERT(curves.GetNumVertices() > 0);
    YBI_ASSERT(curves.GetNumCurves() > 0);

    const Array<Vec3> &positions = curves.GetVertices();
    const Array<float> &widths = curves.GetWidths();
    const uint32_t numVertices = static_cast<uint32_t>(curves.GetNumVertices());

    uint32_t totalSegments = 0;
    for (uint32_t curveIndex = 0; curveIndex < curves.GetNumCurves(); curveIndex++)
    {
        uint32_t start = 0;
        uint32_t count = 0;
        curves.GetCurveRange(curveIndex, start, count);
        if (count >= 4)
        {
            totalSegments += count - 3;
        }
    }

    if (totalSegments == 0)
    {
        return nullptr;
    }

    RTCScene leafScene = rtcNewScene(embreeDevice);
    YBI_ASSERT(leafScene);

    RTCGeometry curveGeometry = rtcNewGeometry(embreeDevice, RTC_GEOMETRY_TYPE_ROUND_BSPLINE_CURVE);
    YBI_ASSERT(curveGeometry);

    Vec4 *vertexBuffer = static_cast<Vec4 *>(rtcSetNewGeometryBuffer(curveGeometry,
                                                                          RTC_BUFFER_TYPE_VERTEX,
                                                                          0,
                                                                          RTC_FORMAT_FLOAT4,
                                                                          sizeof(Vec4),
                                                                          numVertices));
    uint32_t *indexBuffer = static_cast<uint32_t *>(rtcSetNewGeometryBuffer(curveGeometry,
                                                                             RTC_BUFFER_TYPE_INDEX,
                                                                             0,
                                                                             RTC_FORMAT_UINT,
                                                                             sizeof(uint32_t),
                                                                             totalSegments));
    YBI_ASSERT(vertexBuffer);
    YBI_ASSERT(indexBuffer);

    for (uint32_t curveIndex = 0; curveIndex < curves.GetNumCurves(); curveIndex++)
    {
        uint32_t start = 0;
        uint32_t count = 0;
        curves.GetCurveRange(curveIndex, start, count);
        for (uint32_t localVertex = 0; localVertex < count; localVertex++)
        {
            const uint32_t vertexIndex = start + localVertex;
            const Vec3 p = positions[vertexIndex];
            const float width = std::max(0.0f, ResolveCurveWidth(curves, widths, curveIndex, vertexIndex));
            vertexBuffer[vertexIndex] = Vec4(p.x, p.y, p.z, width * 0.5f);
        }
    }

    uint32_t segmentOut = 0;
    for (uint32_t curveIndex = 0; curveIndex < curves.GetNumCurves(); curveIndex++)
    {
        uint32_t start = 0;
        uint32_t count = 0;
        curves.GetCurveRange(curveIndex, start, count);
        if (count < 4)
        {
            continue;
        }
        const uint32_t numSegments = count - 3;
        for (uint32_t segmentIndex = 0; segmentIndex < numSegments; segmentIndex++)
        {
            indexBuffer[segmentOut++] = start + segmentIndex;
        }
    }
    YBI_ASSERT(segmentOut == totalSegments);

    rtcSetGeometryTessellationRate(curveGeometry, 1.0f);
    rtcCommitGeometry(curveGeometry);
    rtcAttachGeometry(leafScene, curveGeometry);
    rtcReleaseGeometry(curveGeometry);
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
        rtcSetGeometryUserData(
            instance,
            reinterpret_cast<void *>(static_cast<uintptr_t>(mesh.refIndex)));
        float transform[12] = {};
        CopyTransformMatrix(mesh.parentFromLocal, transform);
        rtcSetGeometryTransform(instance, 0, RTC_FORMAT_FLOAT3X4_ROW_MAJOR, transform);
        rtcCommitGeometry(instance);
        rtcAttachGeometryByID(embreeScene, instance, mesh.refIndex);
        rtcReleaseGeometry(instance);

        rtcReleaseScene(leafScene);
    }

    for (const Curves &curves : scene->curves)
    {
        if (curves.GetNumVertices() == 0 || curves.GetNumCurves() == 0)
        {
            continue;
        }

        RTCScene leafScene = BuildCurveLeafScene(cpuDevice->embreeDevice, curves);
        if (!leafScene)
        {
            static bool warnedEmptySegments = false;
            if (!warnedEmptySegments)
            {
                fprintf(stderr, "Embree CPU path: curve primitive has no valid segments; skipping.\n");
                warnedEmptySegments = true;
            }
            continue;
        }

        RTCGeometry instance = rtcNewGeometry(cpuDevice->embreeDevice, RTC_GEOMETRY_TYPE_INSTANCE);
        YBI_ASSERT(instance);
        rtcSetGeometryInstancedScene(instance, leafScene);
        rtcSetGeometryUserData(instance, reinterpret_cast<void *>(static_cast<uintptr_t>(UINT32_MAX)));
        float transform[12] = {};
        CopyTransformMatrix(curves.parentFromLocal, transform);
        rtcSetGeometryTransform(instance, 0, RTC_FORMAT_FLOAT3X4_ROW_MAJOR, transform);
        rtcCommitGeometry(instance);
        rtcAttachGeometry(embreeScene, instance);
        rtcReleaseGeometry(instance);

        rtcReleaseScene(leafScene);
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
        uint32_t childRefIndex = UINT32_MAX;
        if (!childScene->meshes.empty())
        {
            childRefIndex = childScene->meshes[0].refIndex;
        }
        rtcSetGeometryUserData(
            childInstance,
            reinterpret_cast<void *>(static_cast<uintptr_t>(childRefIndex)));
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

} // namespace ybi

#endif
