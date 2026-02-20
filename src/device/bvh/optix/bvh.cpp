#include "device/cuda_device.h"

#include "scene/scene.h"

#include <cstring>
#include <vector>

#if (OPTIX_VERSION >= 90000)
#include "device/bvh/optix/clusters.cuh"
#endif

YBI_NAMESPACE_BEGIN

#if defined(WITH_CUDA) && defined(WITH_OPTIX)

static OptixAccelBuildOptions GetDefaultBuildOptions(uint32_t numMotionKeys)
{
    OptixAccelBuildOptions options = {};
    options.buildFlags = OPTIX_BUILD_FLAG_ALLOW_COMPACTION | OPTIX_BUILD_FLAG_PREFER_FAST_TRACE;
    options.operation = OPTIX_BUILD_OPERATION_BUILD;
    options.motionOptions.numKeys = numMotionKeys;
    options.motionOptions.flags = 0;
    options.motionOptions.timeBegin = 0.f;
    options.motionOptions.timeEnd = 1.f;
    return options;
}

static bool IsIdentityTransform(const float3x4 &transform)
{
    const float identity[3][4] = {
        {1.0f, 0.0f, 0.0f, 0.0f},
        {0.0f, 1.0f, 0.0f, 0.0f},
        {0.0f, 0.0f, 1.0f, 0.0f},
    };
    return memcmp(transform.m, identity, sizeof(identity)) == 0;
}

static void CopyTransformMatrix(const float3x4 &transform, float out[12])
{
    memcpy(out, transform.m, sizeof(float) * 12);
}

static void SetOptixInstanceDefaults(OptixInstance &instance)
{
    memset(&instance, 0, sizeof(instance));
    instance.visibilityMask = 0xFFu;
    instance.flags = OPTIX_INSTANCE_FLAG_NONE;
}

static OptixTraversableHandle BuildOptixBVH(CUDADevice *cudaDevice,
                                            CUDAMemoryArena &deviceArena,
                                            OptixAccelBuildOptions buildOptions,
                                            const OptixBuildInput *buildInputs,
                                            uint32_t numBuildInputs)
{
    YBI_ASSERT(buildInputs);
    YBI_ASSERT(numBuildInputs > 0);

    OptixAccelBufferSizes sizes = {};
    OPTIX_ASSERT(optixAccelComputeMemoryUsage(
        cudaDevice->optixDeviceContext, &buildOptions, buildInputs, numBuildInputs, &sizes));

    cudaDevice->bvhTotalAllocated += sizes.outputSizeInBytes;

    const bool useCompaction = (buildOptions.buildFlags & OPTIX_BUILD_FLAG_ALLOW_COMPACTION) != 0;
    DeviceMemoryView<uint8_t> tempBuffer = deviceArena.PushArray<uint8_t>(
        util::AlignUp(sizes.tempSizeInBytes + sizeof(uint64_t), sizeof(uint64_t)),
        OPTIX_ACCEL_BUFFER_BYTE_ALIGNMENT);

    DeviceMemoryView<uint8_t> outputBuffer =
        useCompaction ? deviceArena.PushArray<uint8_t>(sizes.outputSizeInBytes,
                                                       OPTIX_ACCEL_BUFFER_BYTE_ALIGNMENT)
                      : cudaDevice->Alloc<uint8_t>(sizes.outputSizeInBytes);

    uint64_t compactedSizeAddress =
        util::AlignUp(CUdeviceptr(tempBuffer.data()) + sizes.tempSizeInBytes, sizeof(uint64_t));
    OptixAccelEmitDesc emittedProperties = {};
    emittedProperties.type = OPTIX_PROPERTY_TYPE_COMPACTED_SIZE;
    emittedProperties.result = CUdeviceptr(compactedSizeAddress);

    OptixTraversableHandle outputHandle = {};
    OPTIX_ASSERT(optixAccelBuild(cudaDevice->optixDeviceContext,
                                 0,
                                 &buildOptions,
                                 buildInputs,
                                 numBuildInputs,
                                 CUdeviceptr(tempBuffer.data()),
                                 sizes.tempSizeInBytes,
                                 CUdeviceptr(outputBuffer.data()),
                                 sizes.outputSizeInBytes,
                                 &outputHandle,
                                 useCompaction ? &emittedProperties : nullptr,
                                 useCompaction ? 1 : 0));

    CUDA_ASSERT(cuStreamSynchronize(0));

    if (useCompaction)
    {
        uint64_t compactedSize = sizes.outputSizeInBytes;
        CUDA_ASSERT(
            cuMemcpyDtoH(&compactedSize, CUdeviceptr(compactedSizeAddress), sizeof(uint64_t)));

        DeviceMemoryView<uint8_t> compactedOutputBuffer =
            cudaDevice->Alloc<uint8_t>(compactedSize);

        OptixTraversableHandle tempOutputHandle = {};
        OPTIX_ASSERT(optixAccelCompact(cudaDevice->optixDeviceContext,
                                       0,
                                       outputHandle,
                                       CUdeviceptr(compactedOutputBuffer.data()),
                                       compactedSize,
                                       &tempOutputHandle));
        CUDA_ASSERT(cuStreamSynchronize(0));
        outputHandle = tempOutputHandle;

        cudaDevice->bvhTotalAllocated -= sizes.outputSizeInBytes;
        cudaDevice->bvhTotalAllocated += compactedSize;
    }

    return outputHandle;
}

static OptixBuildInput GetOptiXTriangleBuildInput(CUDADevice *cudaDevice,
                                                  HostMemoryArena &hostArena,
                                                  CUDAMemoryArena &deviceArena,
                                                  Mesh &mesh,
                                                  uint32_t numMotionKeys,
                                                  OptixAccelBuildOptions &options)
{
    (void)cudaDevice;
    (void)options;

    const uint32_t numVertices = mesh.positions.size();
    const uint32_t numIndices = mesh.indices.size();

    MemoryView<CUdeviceptr> vertexBuffers = hostArena.PushArray<CUdeviceptr>(numMotionKeys);
    MemoryView<float3> hostVertices = hostArena.PushArray<float3>(numVertices * numMotionKeys);

    DeviceMemoryView<float3> deviceVertices =
        deviceArena.PushArray<float3>(numVertices * numMotionKeys);
    DeviceMemoryView<int> deviceIndices = deviceArena.PushArray<int>(numIndices);
    DeviceMemoryView<float> preTransform = deviceArena.PushArray<float>(12);

    for (uint32_t step = 0; step < numMotionKeys; step++)
    {
        CUdeviceptr dst = (CUdeviceptr)(deviceVertices.data() + step * numVertices);
        vertexBuffers[step] = dst;

        memcpy(hostVertices.data() + step * numVertices,
               mesh.positions.data() + step * numVertices,
               sizeof(float3) * numVertices);
    }

    float hostPreTransform[12] = {};
    CopyTransformMatrix(mesh.parentFromLocal, hostPreTransform);

    CUDA_ASSERT(cuMemcpyHtoD(
        CUdeviceptr(deviceVertices.data()), hostVertices.data(), deviceVertices.numBytes()));
    CUDA_ASSERT(cuMemcpyHtoD(
        CUdeviceptr(deviceIndices.data()), mesh.indices.data(), deviceIndices.numBytes()));
    CUDA_ASSERT(cuMemcpyHtoD(
        CUdeviceptr(preTransform.data()), hostPreTransform, sizeof(hostPreTransform)));

    MemoryView<unsigned int> geometryFlags = hostArena.PushArray<unsigned int>(1);
    geometryFlags[0] = OPTIX_GEOMETRY_FLAG_DISABLE_ANYHIT;

    OptixBuildInput input = {};
    input.type = OPTIX_BUILD_INPUT_TYPE_TRIANGLES;
    OptixBuildInputTriangleArray &triangleArray = input.triangleArray;
    triangleArray = {};
    triangleArray.vertexBuffers = vertexBuffers.data();
    triangleArray.numVertices = numVertices;
    triangleArray.vertexFormat = OPTIX_VERTEX_FORMAT_FLOAT3;
    triangleArray.vertexStrideInBytes = sizeof(float3);
    triangleArray.indexBuffer = CUdeviceptr(deviceIndices.data());
    triangleArray.numIndexTriplets = numIndices / 3;
    triangleArray.indexFormat = OPTIX_INDICES_FORMAT_UNSIGNED_INT3;
    triangleArray.indexStrideInBytes = sizeof(int) * 3;
    triangleArray.preTransform = CUdeviceptr(preTransform.data());
    triangleArray.transformFormat = OPTIX_TRANSFORM_FORMAT_MATRIX_FLOAT12;
    triangleArray.flags = geometryFlags.data();
    triangleArray.numSbtRecords = 1;

    return input;
}

static OptixBuildInput GetOptiXCurveBuildInput(CUDADevice *cudaDevice,
                                               HostMemoryArena &hostArena,
                                               CUDAMemoryArena &deviceArena,
                                               Curves &curves,
                                               uint32_t numMotionKeys,
                                               OptixAccelBuildOptions &options)
{
    (void)cudaDevice;
    (void)options;

    const uint32_t numVertices = (uint32_t)curves.GetNumVertices();
    const uint32_t totalNumSegments = (uint32_t)curves.GetNumSegments();
    const uint32_t totalNumVertices = numVertices * numMotionKeys;

    MemoryView<CUdeviceptr> vertexBuffers = hostArena.PushArray<CUdeviceptr>(numMotionKeys);
    MemoryView<CUdeviceptr> widthBuffers = hostArena.PushArray<CUdeviceptr>(numMotionKeys);
    MemoryView<uint32_t> indexBuffer = hostArena.PushArray<uint32_t>(totalNumSegments);
    MemoryView<float3> hostVertices = hostArena.PushArray<float3>(totalNumVertices);
    MemoryView<float> hostWidths = hostArena.PushArray<float>(totalNumVertices);

    DeviceMemoryView<float3> deviceVertices = deviceArena.PushArray<float3>(totalNumVertices);
    DeviceMemoryView<int> deviceIndices = deviceArena.PushArray<int>(totalNumSegments);
    DeviceMemoryView<float> deviceWidths = deviceArena.PushArray<float>(totalNumVertices);

    size_t segmentIndexOut = 0;
    for (size_t curveIndex = 0; curveIndex < curves.GetNumCurves(); curveIndex++)
    {
        const int segmentStart = curves.GetCurveKeyStart(curveIndex);
        const int numSegments = curves.GetCurveNumSegments(curveIndex);
        for (int segmentIndex = segmentStart; segmentIndex < segmentStart + numSegments;
             segmentIndex++, segmentIndexOut++)
        {
            indexBuffer[segmentIndexOut] = (uint32_t)segmentIndex;
        }
    }

    for (uint32_t step = 0; step < numMotionKeys; step++)
    {
        CUdeviceptr dst = (CUdeviceptr)(deviceVertices.data() + step * numVertices);
        vertexBuffers[step] = dst;

        dst = (CUdeviceptr)(deviceWidths.data() + step * numVertices);
        widthBuffers[step] = dst;

        const Array<float3> &positions = curves.GetVertices();
        const Array<float> &widths = curves.GetWidths();
        memcpy(hostVertices.data() + step * numVertices,
               positions.data(),
               sizeof(float3) * numVertices);
        memcpy(hostWidths.data() + step * numVertices, widths.data(), sizeof(float) * numVertices);
    }

    CUDA_ASSERT(cuMemcpyHtoD(
        CUdeviceptr(deviceVertices.data()), hostVertices.data(), deviceVertices.numBytes()));
    CUDA_ASSERT(cuMemcpyHtoD(
        CUdeviceptr(deviceIndices.data()), indexBuffer.data(), deviceIndices.numBytes()));
    CUDA_ASSERT(cuMemcpyHtoD(
        CUdeviceptr(deviceWidths.data()), hostWidths.data(), deviceWidths.numBytes()));

    if (!IsIdentityTransform(curves.parentFromLocal))
    {
        // TODO: OptixBuildInputCurveArray has no preTransform; apply transforms through IAS
        // instances.
    }

    OptixBuildInput input = {};
    input.type = OPTIX_BUILD_INPUT_TYPE_CURVES;
    OptixBuildInputCurveArray &curveArray = input.curveArray;
    curveArray = {};
    curveArray.curveType = OPTIX_PRIMITIVE_TYPE_ROUND_CUBIC_BSPLINE;
    curveArray.numPrimitives = totalNumSegments;
    curveArray.vertexBuffers = vertexBuffers.data();
    curveArray.numVertices = numVertices;
    curveArray.vertexStrideInBytes = sizeof(float3);
    curveArray.widthBuffers = widthBuffers.data();
    curveArray.widthStrideInBytes = sizeof(float);
    curveArray.normalBuffers = nullptr;
    curveArray.normalStrideInBytes = 0;
    curveArray.indexBuffer = CUdeviceptr(deviceIndices.data());
    curveArray.indexStrideInBytes = sizeof(int);
    curveArray.flag = OPTIX_GEOMETRY_FLAG_NONE;
    curveArray.primitiveIndexOffset = 0;
    curveArray.endcapFlags = OPTIX_CURVE_ENDCAP_DEFAULT;

    return input;
}

static OptixBuildInput GetOptiXInstanceBuildInput(CUDAMemoryArena &deviceArena,
                                                  MemoryView<OptixInstance> hostInstances)
{
    YBI_ASSERT(hostInstances.size() > 0);

    DeviceMemoryView<OptixInstance> deviceInstances =
        deviceArena.PushArray<OptixInstance>(hostInstances.size());
    CUDA_ASSERT(cuMemcpyHtoD(
        CUdeviceptr(deviceInstances.data()), hostInstances.data(), deviceInstances.numBytes()));

    OptixBuildInput input = {};
    input.type = OPTIX_BUILD_INPUT_TYPE_INSTANCES;
    input.instanceArray.instances = CUdeviceptr(deviceInstances.data());
    input.instanceArray.numInstances = (unsigned int)hostInstances.size();

    return input;
}

static OptixTraversableHandle
BuildSceneGeometryGAS(CUDADevice *cudaDevice, HostMemoryArena &hostArena, Scene *scene)
{
    const uint32_t numMotionKeys = 1;
    OptixAccelBuildOptions options = GetDefaultBuildOptions(numMotionKeys);

    std::vector<OptixBuildInput> buildInputs;
    buildInputs.reserve(scene->meshes.size() + scene->curves.size());

    for (Mesh &mesh : scene->meshes)
    {
        buildInputs.push_back(GetOptiXTriangleBuildInput(
            cudaDevice, hostArena, *cudaDevice->deviceArena, mesh, numMotionKeys, options));
    }

    for (Curves &curves : scene->curves)
    {
        buildInputs.push_back(GetOptiXCurveBuildInput(
            cudaDevice, hostArena, *cudaDevice->deviceArena, curves, numMotionKeys, options));
    }

    if (buildInputs.empty())
    {
        return {};
    }

    return BuildOptixBVH(cudaDevice,
                         *cudaDevice->deviceArena,
                         options,
                         buildInputs.data(),
                         (uint32_t)buildInputs.size());
}

static OptixTraversableHandle BuildSceneInstanceIAS(CUDADevice *cudaDevice,
                                                    HostMemoryArena &hostArena,
                                                    Scene *scene,
                                                    OptixTraversableHandle localGeometryHandle)
{
    const bool hasLocalGeometry = localGeometryHandle != OptixTraversableHandle(0);
    const size_t numInstances = scene->instances.size() + (hasLocalGeometry ? 1 : 0);
    if (numInstances == 0)
    {
        return {};
    }

    MemoryView<OptixInstance> hostInstances = hostArena.PushArray<OptixInstance>(numInstances);

    size_t instanceWriteIndex = 0;
    if (hasLocalGeometry)
    {
        OptixInstance &localGeometryInstance = hostInstances[instanceWriteIndex++];
        SetOptixInstanceDefaults(localGeometryInstance);

        static const float3x4 identity = {
            1.0f,
            0.0f,
            0.0f,
            0.0f,
            0.0f,
            1.0f,
            0.0f,
            0.0f,
            0.0f,
            0.0f,
            1.0f,
            0.0f,
        };
        CopyTransformMatrix(identity, localGeometryInstance.transform);
        localGeometryInstance.instanceId = (unsigned int)(instanceWriteIndex - 1);
        localGeometryInstance.sbtOffset = 0;
        localGeometryInstance.traversableHandle = localGeometryHandle;
    }

    for (const Instance &instance : scene->instances)
    {
        YBI_ASSERT(instance.childSceneIndex < scene->childScenes.size());
        Scene *childScene = scene->childScenes[instance.childSceneIndex];
        YBI_ASSERT(childScene);
        YBI_ASSERT(childScene->bvhHandle != 0);

        OptixInstance &optixInstance = hostInstances[instanceWriteIndex++];
        SetOptixInstanceDefaults(optixInstance);
        CopyTransformMatrix(instance.parentFromLocal, optixInstance.transform);
        optixInstance.instanceId = (unsigned int)(instanceWriteIndex - 1);
        optixInstance.sbtOffset = 0;
        optixInstance.traversableHandle = (OptixTraversableHandle)childScene->bvhHandle;
    }

    OptixBuildInput buildInput =
        GetOptiXInstanceBuildInput(*cudaDevice->deviceArena, hostInstances);
    const uint32_t numMotionKeys = 1;
    OptixAccelBuildOptions options = GetDefaultBuildOptions(numMotionKeys);

    return BuildOptixBVH(cudaDevice, *cudaDevice->deviceArena, options, &buildInput, 1);
}

OptixTraversableHandle
BuildTriangleGASFromMesh(CUDADevice *cudaDevice, HostMemoryArena &hostArena, Mesh &mesh)
{
    CUDA_ASSERT(cuCtxPushCurrent(cudaDevice->cudaContext));

    const uint32_t numMotionKeys = 1;
    OptixAccelBuildOptions options = GetDefaultBuildOptions(numMotionKeys);

    OptixBuildInput buildInput = GetOptiXTriangleBuildInput(
        cudaDevice, hostArena, *cudaDevice->deviceArena, mesh, numMotionKeys, options);
    OptixTraversableHandle handle =
        BuildOptixBVH(cudaDevice, *cudaDevice->deviceArena, options, &buildInput, 1);

    CUDA_ASSERT(cuCtxPopCurrent(0));
    return handle;
}

OptixTraversableHandle
BuildClusterGASFromMesh(CUDADevice *cudaDevice, HostMemoryArena &hostArena, const Mesh &mesh)
{
#if (OPTIX_VERSION >= 90000)
    CUDA_ASSERT(cuCtxPushCurrent(cudaDevice->cudaContext));
    MeshletClustersResult clusterResult = ClusterizeTest(mesh);
    YBI_ASSERT(clusterResult.meshletCount > 0);
    size_t totalOutputSize = 0;
    OptixTraversableHandle handle = {};
    BuildMeshCLAS(cudaDevice, hostArena, mesh, clusterResult, totalOutputSize, handle);
    (void)totalOutputSize;
    CUDA_ASSERT(cuCtxPopCurrent(0));
    return handle;
#else
    (void)cudaDevice;
    (void)hostArena;
    (void)mesh;
    return {};
#endif
}

OptixTraversableHandle
BuildCurveGASFromCurves(CUDADevice *cudaDevice, HostMemoryArena &hostArena, Curves &curves)
{
    CUDA_ASSERT(cuCtxPushCurrent(cudaDevice->cudaContext));

    const uint32_t numMotionKeys = 1;
    OptixAccelBuildOptions options = GetDefaultBuildOptions(numMotionKeys);

    OptixBuildInput buildInput = GetOptiXCurveBuildInput(
        cudaDevice, hostArena, *cudaDevice->deviceArena, curves, numMotionKeys, options);
    OptixTraversableHandle handle =
        BuildOptixBVH(cudaDevice, *cudaDevice->deviceArena, options, &buildInput, 1);

    CUDA_ASSERT(cuCtxPopCurrent(0));
    return handle;
}

// TODO: handle the case where cuda is enabled but optix isn't
void CUDADevice::BuildBVH(ScenePool *scenePool)
{
    YBI_ASSERT(scenePool);

    CUDA_ASSERT(cuCtxPushCurrent(cudaContext));

    HostMemoryArena hostArena;

    for (const std::unique_ptr<Scene> &scenePtr : scenePool->scenes)
    {
        Scene *scene = scenePtr.get();
        YBI_ASSERT(scene);

        scene->bvhHandle = 0;

        OptixTraversableHandle localGeometryHandle = BuildSceneGeometryGAS(this, hostArena, scene);
        OptixTraversableHandle sceneHandle = localGeometryHandle;

        if (!scene->instances.empty())
        {
            sceneHandle = BuildSceneInstanceIAS(this, hostArena, scene, localGeometryHandle);
        }

        scene->bvhHandle = (Scene::BVHHandle)sceneHandle;

        hostArena.Clear();
        deviceArena->Clear();
    }

    CUDA_ASSERT(cuCtxPopCurrent(0));
}

#endif

YBI_NAMESPACE_END
