#pragma once

#include "cuda.h"
#include "cuda_runtime_api.h"
#include "device/cuda_assert.h"
#include "device/cuda_device_memory_arena.h"
#include "optix_types.h"
#include "scene/micropolygon_mesh.h"
#include "scene/scene.h"
#include "util/array.h"
#include "util/assert.h"
#include "util/base.h"
#include "util/float3.h"
#include "util/host_memory_arena.h"

#include <cuda_runtime.h>
#include <optix_function_table_definition.h>
#include <optix_stubs.h>
#include <string>

YBI_NAMESPACE_BEGIN

#ifdef WITH_CUDA

#ifdef WITH_OPTIX

struct ClusterAccelerationStructureLimits
{
    uint32_t maxTrianglesPerCluster;
    uint32_t maxVerticesPerCluster;
};

static OptixDeviceContext InitializeOptix(CUcontext cudaContext)
{
    OPTIX_ASSERT(optixInit());
    OptixDeviceContextOptions contextOptions = {};
    contextOptions.logCallbackFunction =
        [](unsigned int level, const char *tag, const char *message, void *cbdata) {
            std::string type = {};

            switch (level)
            {
                case 1:
                    type = "Fatal Error";
                    break;
                case 2:
                    type = "Error";
                    break;
                case 3:
                    type = "Warning";
                    break;
                case 4:
                    type = "Status";
                    break;
                default:
                    break;
            }

            // Print("Optix %S: %s\n", type, message);
        };
    contextOptions.logCallbackLevel = 4;
    contextOptions.validationMode = OPTIX_DEVICE_CONTEXT_VALIDATION_MODE_ALL;

    OptixDeviceContext optixDeviceContext;
    OPTIX_ASSERT(optixDeviceContextCreate(cudaContext, &contextOptions, &optixDeviceContext));
    OPTIX_ASSERT(optixDeviceContextSetLogCallback(optixDeviceContext,
                                                  contextOptions.logCallbackFunction,
                                                  contextOptions.logCallbackData,
                                                  contextOptions.logCallbackLevel));
    return optixDeviceContext;
}

// TODO: gpu kernels shouldn't be in the same file as host code
extern "C" __global__ void
Test(uint8_t *indexBuffer, float3 *vertexBuffer, ClusterAccelerationStructureLimits limits)
{
    uint32_t clusterId = 0;
    uint32_t triangleCount = 0;
    uint32_t vertexCount = 0;
    YBI_ASSERT(triangleCount <= limits.maxTrianglesPerCluster);
    YBI_ASSERT(vertexCount <= limits.maxVerticesPerCluster);

    OptixClusterAccelBuildInputTrianglesArgs args = {};
    args.clusterId = clusterId;
    args.clusterFlags = 0;
    args.triangleCount = triangleCount;
    args.vertexCount = vertexCount;
    args.positionTruncateBitCount = 0;
    args.indexFormat = OPTIX_CLUSTER_ACCEL_INDICES_FORMAT_8BIT;
    args.indexBuffer = (CUdeviceptr)indexBuffer;
    args.vertexBuffer = (CUdeviceptr)vertexBuffer;
}

struct CUDADevice
{
    CUcontext cudaContext;
    OptixDeviceContext optixDeviceContext;
    size_t totalAllocated;
    size_t bvhTotalAllocated;

#if (OPTIX_VERSION >= 90000)
    DeviceMemoryView<OptixTraversableHandle> gridClusterTemplateOutputHandles;
#endif

    CUDADevice();

    template <typename T>
    DeviceMemoryView<T> Alloc(size_t size);
    template <typename T>
    void Free(DeviceMemoryView<T> &view);
    bool SupportsGrids() const;

    void CreateGridClusterTemplates();
};

CUDADevice::CUDADevice() : totalAllocated(0), bvhTotalAllocated()
{
    cuInit(0);
    CUdevice device;
    cuDeviceGet(&device, 0);
    cuDevicePrimaryCtxRetain(&cudaContext, device);

    optixDeviceContext = InitializeOptix(cudaContext);
}

template <typename T>
DeviceMemoryView<T> CUDADevice::Alloc(size_t count)
{
    YBI_ASSERT(count != 0);
    size_t size = sizeof(T) * count;
    totalAllocated += size;
    CUdeviceptr ptr;
    CUDA_ASSERT(cuMemAlloc(&ptr, size));

    return {(T *)ptr, count};
}

template <typename T>
void CUDADevice::Free(DeviceMemoryView<T> &view)
{
    if (view.data() == nullptr || view.size() == 0)
        return;
    totalAllocated -= view.numBytes();
    CUDA_ASSERT(cuMemFree(CUdeviceptr(view.data())));
    view = {};
}

bool CUDADevice::SupportsGrids() const
{
#if (OPTIX_VERSION >= 90000)
    return true;
#else
    return false;
#endif
}

void CUDADevice::CreateGridClusterTemplates()
{
#if (OPTIX_VERSION >= 90000)
    const int minDim = 2;
    const int maxDim = 8;
    const int widthDim = maxDim - minDim + 1;
    const int numGrids = widthDim * widthDim;

    unsigned int maxEdges;
    OPTIX_ASSERT(
        optixDeviceContextGetProperty(optixDeviceContext,
                                      OPTIX_DEVICE_PROPERTY_LIMIT_MAX_STRUCTURED_GRID_RESOLUTION,
                                      &maxEdges,
                                      sizeof(unsigned int)));
    YBI_ASSERT(maxDim <= (int)maxEdges);

    OptixClusterAccelBuildInput input = {};
    input.type = OPTIX_CLUSTER_ACCEL_BUILD_TYPE_TEMPLATES_FROM_GRIDS;
    input.grids.flags = OPTIX_CLUSTER_ACCEL_BUILD_FLAG_PREFER_FAST_TRACE;
    input.grids.maxArgCount = (unsigned int)numGrids;
    input.grids.vertexFormat = OPTIX_VERTEX_FORMAT_FLOAT3;
    input.grids.maxSbtIndexValue = 0;
    input.grids.maxWidth = (unsigned int)maxDim;
    input.grids.maxHeight = (unsigned int)maxDim;

    OptixClusterAccelBuildInputGridsArgs hostGridArgs[49];
    int idx = 0;
    for (int w = minDim; w <= maxDim; w++)
    {
        for (int h = minDim; h <= maxDim; h++)
        {
            OptixClusterAccelBuildInputGridsArgs &arg = hostGridArgs[idx];
            arg = {};
            arg.baseClusterId = 0;
            arg.clusterFlags = 0;
            arg.positionTruncateBitCount = 0;
            arg.dimensions[0] = (unsigned char)w;
            arg.dimensions[1] = (unsigned char)h;
            idx++;
        }
    }

    DeviceMemoryView<OptixClusterAccelBuildInputGridsArgs> args =
        Alloc<OptixClusterAccelBuildInputGridsArgs>((size_t)numGrids);
    DeviceMemoryView<uint32_t> argsCount = Alloc<uint32_t>(1);
    const uint32_t argsCountVal = (uint32_t)numGrids;
    CUDA_ASSERT(cuMemcpyHtoD(CUdeviceptr(args.data()),
                             hostGridArgs,
                             (size_t)numGrids * sizeof(OptixClusterAccelBuildInputGridsArgs)));
    CUDA_ASSERT(cuMemcpyHtoD(CUdeviceptr(argsCount.data()), &argsCountVal, sizeof(uint32_t)));

    OptixAccelBufferSizes bufferSizes = {};
    OPTIX_ASSERT(
        optixClusterAccelComputeMemoryUsage(optixDeviceContext,
                                            OPTIX_CLUSTER_ACCEL_BUILD_MODE_IMPLICIT_DESTINATIONS,
                                            &input,
                                            &bufferSizes));

    YBI_ASSERT(bufferSizes.tempSizeInBytes);
    DeviceMemoryView<char> temp = Alloc<char>(bufferSizes.tempSizeInBytes);

    YBI_ASSERT(bufferSizes.outputSizeInBytes);
    DeviceMemoryView<char> output = Alloc<char>(bufferSizes.outputSizeInBytes);

    gridClusterTemplateOutputHandles = Alloc<OptixTraversableHandle>((size_t)numGrids);

    OptixClusterAccelBuildModeDesc desc = {};
    desc.mode = OPTIX_CLUSTER_ACCEL_BUILD_MODE_IMPLICIT_DESTINATIONS;
    desc.implicitDest.tempBuffer = CUdeviceptr(temp.data());
    desc.implicitDest.tempBufferSizeInBytes = bufferSizes.tempSizeInBytes;
    desc.implicitDest.outputBuffer = CUdeviceptr(output.data());
    desc.implicitDest.outputBufferSizeInBytes = bufferSizes.outputSizeInBytes;
    desc.implicitDest.outputHandlesBuffer = CUdeviceptr(gridClusterTemplateOutputHandles.data());
    desc.implicitDest.outputHandlesStrideInBytes = sizeof(OptixTraversableHandle);

    OPTIX_ASSERT(optixClusterAccelBuild(optixDeviceContext,
                                        0,
                                        &desc,
                                        &input,
                                        CUdeviceptr(args.data()),
                                        CUdeviceptr(argsCount.data()),
                                        sizeof(OptixClusterAccelBuildInputGridsArgs)));

    CUDA_ASSERT(cuStreamSynchronize(0));
    Free(temp);
#endif
}

static OptixTraversableHandle BuildOptixBVH(CUDADevice *cudaDevice,
                                            CUDAMemoryArena &deviceArena,
                                            OptixAccelBuildOptions buildOptions,
                                            OptixBuildInput buildInput)
{
    OptixAccelBufferSizes sizes = {};
    OPTIX_ASSERT(optixAccelComputeMemoryUsage(
        cudaDevice->optixDeviceContext, &buildOptions, &buildInput, 1, &sizes));

    cudaDevice->bvhTotalAllocated += sizes.outputSizeInBytes;

    // Build BVH
    bool useCompaction = buildOptions.buildFlags & OPTIX_BUILD_FLAG_ALLOW_COMPACTION;
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
    OptixTraversableHandle outputHandle;
    OPTIX_ASSERT(optixAccelBuild(cudaDevice->optixDeviceContext,
                                 0,
                                 &buildOptions,
                                 &buildInput,
                                 1,
                                 CUdeviceptr(tempBuffer.data()),
                                 sizes.tempSizeInBytes,
                                 CUdeviceptr(outputBuffer.data()),
                                 sizes.outputSizeInBytes,
                                 &outputHandle,
                                 useCompaction ? &emittedProperties : 0,
                                 useCompaction ? 1 : 0));

    CUDA_ASSERT(cuStreamSynchronize(0));

    if (useCompaction)
    {
        uint64_t compactedSize = sizes.outputSizeInBytes;
        CUDA_ASSERT(
            cuMemcpyDtoH(&compactedSize, CUdeviceptr(compactedSizeAddress), sizeof(uint64_t)));

        if (compactedSize < sizes.outputSizeInBytes)
        {
            DeviceMemoryView<uint8_t> compactedOutputBuffer =
                cudaDevice->Alloc<uint8_t>(compactedSize);

            OptixTraversableHandle tempOutputHandle;
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
    const uint32_t numVertices = mesh.positions.size();
    const uint32_t numIndices = mesh.indices.size();

    MemoryView<CUdeviceptr> vertexBuffers = hostArena.PushArray<CUdeviceptr>(numMotionKeys);
    MemoryView<float3> hostVertices = hostArena.PushArray<float3>(numVertices * numMotionKeys);

    DeviceMemoryView<float3> deviceVertices =
        deviceArena.PushArray<float3>(numVertices * numMotionKeys);
    DeviceMemoryView<int> deviceIndices = deviceArena.PushArray<int>(numIndices);

    for (uint32_t step = 0; step < numMotionKeys; step++)
    {
        CUdeviceptr dst = (CUdeviceptr)((deviceVertices + step * numVertices).data());
        vertexBuffers[step] = dst;

        util::Copy(
            hostVertices + step * numVertices, mesh.positions + step * numVertices, numVertices);
    }
    CUDA_ASSERT(cuMemcpyHtoD(
        CUdeviceptr(deviceVertices.data()), hostVertices.data(), deviceVertices.numBytes()));
    CUDA_ASSERT(cuMemcpyHtoD(
        CUdeviceptr(deviceIndices.data()), mesh.indices.data(), deviceIndices.numBytes()));

    unsigned int flags = OPTIX_GEOMETRY_FLAG_REQUIRE_SINGLE_ANYHIT_CALL;
    OptixBuildInput input = {};
    input.type = OPTIX_BUILD_INPUT_TYPE_TRIANGLES;
    OptixBuildInputTriangleArray &triangleArray = input.triangleArray;
    triangleArray = {};
    triangleArray.vertexBuffers = vertexBuffers.data();
    triangleArray.numVertices = numVertices;
    triangleArray.vertexFormat = OPTIX_VERTEX_FORMAT_FLOAT3;
    triangleArray.vertexStrideInBytes = 0;
    triangleArray.indexBuffer = CUdeviceptr(deviceIndices.data());
    triangleArray.numIndexTriplets = numIndices / 3;
    triangleArray.indexFormat = OPTIX_INDICES_FORMAT_UNSIGNED_INT3;
    triangleArray.indexStrideInBytes = 0;
    triangleArray.flags = &flags;
    triangleArray.numSbtRecords = 1;

    BuildOptixBVH(cudaDevice, deviceArena, options, input);
    return input;
}

// TODO: handle the case where cuda is enabled but optix isn't
void BuildBVH(CUDADevice *cudaDevice, Scene *scene)
{
    CUDA_ASSERT(cuCtxPushCurrent(cudaDevice->cudaContext));

    HostMemoryArena hostArena;
    CUDAMemoryArena deviceArena;

    for (Mesh &mesh : scene->meshes)
    {
        uint32_t numMotionKeys = 1;

        OptixAccelBuildOptions options = {};
        options.buildFlags =
            OPTIX_BUILD_FLAG_ALLOW_COMPACTION | OPTIX_BUILD_FLAG_PREFER_FAST_TRACE;
        options.operation = OPTIX_BUILD_OPERATION_BUILD;
        options.motionOptions.numKeys = numMotionKeys;
        options.motionOptions.flags = 0;
        options.motionOptions.timeBegin = 0.f;
        options.motionOptions.timeEnd = 1.f;

        OptixBuildInput buildInput = GetOptiXTriangleBuildInput(
            cudaDevice, hostArena, deviceArena, mesh, numMotionKeys, options);

        hostArena.Clear();
        deviceArena.Clear();
    }

    for (Curves &curve : scene->curves)
    {
        uint32_t numMotionKeys = 1;

        OptixAccelBuildOptions options = {};
        options.buildFlags =
            OPTIX_BUILD_FLAG_ALLOW_COMPACTION | OPTIX_BUILD_FLAG_PREFER_FAST_TRACE;
        options.operation = OPTIX_BUILD_OPERATION_BUILD;
        options.motionOptions.numKeys = numMotionKeys;
        options.motionOptions.flags = 0;
        options.motionOptions.timeBegin = 0.f;
        options.motionOptions.timeEnd = 1.f;

        size_t totalNumSegments = curve.GetNumSegments();
        MemoryView<CUdeviceptr> vertexBuffers = hostArena.PushArray<CUdeviceptr>(numMotionKeys);
        MemoryView<CUdeviceptr> widthBuffers = hostArena.PushArray<CUdeviceptr>(numMotionKeys);
        MemoryView<uint32_t> indexBuffer = hostArena.PushArray<uint32_t>(totalNumSegments);

        size_t bufferIndex = 0;
        for (size_t curveIndex = 0; curveIndex < curve.GetNumCurves(); curveIndex++)
        {
            int segmentStart = curve.GetCurveKeyStart(curveIndex);
            int numSegments = curve.GetCurveNumSegments(curveIndex);
            for (int segmentIndex = segmentStart; segmentIndex < segmentStart + numSegments;
                 segmentIndex++, bufferIndex++)
            {
                indexBuffer[bufferIndex] = segmentIndex;
            }
        }

        uint32_t numVertices = curve.GetNumVertices();
        uint32_t totalNumVertices = numVertices * numMotionKeys;
        MemoryView<float3> hostVertices = hostArena.PushArray<float3>(totalNumVertices);
        MemoryView<float> hostWidths = hostArena.PushArray<float>(totalNumVertices);

        size_t vertexSize = sizeof(float3) * totalNumVertices;
        size_t indexSize = sizeof(int) * totalNumSegments;

        DeviceMemoryView<float3> deviceVertices = deviceArena.PushArray<float3>(totalNumVertices);
        DeviceMemoryView<int> deviceIndices = deviceArena.PushArray<int>(totalNumSegments);
        DeviceMemoryView<float> deviceWidths = deviceArena.PushArray<float>(totalNumVertices);

        for (uint32_t step = 0; step < numMotionKeys; step++)
        {
            CUdeviceptr dst = (CUdeviceptr)((deviceVertices + step * numVertices).data());
            vertexBuffers[step] = dst;

            dst = (CUdeviceptr)((deviceWidths + step * numVertices).data());
            widthBuffers[step] = dst;

            const Array<float3> &positions = curve.GetVertices();
            const Array<float> &widths = curve.GetWidths();
            util::Copy(hostVertices + step * numVertices, positions, numVertices);
            util::Copy(hostWidths + step * numVertices, widths, numVertices);
        }

        CUDA_ASSERT(cuMemcpyHtoD(
            CUdeviceptr(deviceVertices.data()), hostVertices.data(), deviceVertices.numBytes()));
        CUDA_ASSERT(cuMemcpyHtoD(
            CUdeviceptr(deviceIndices.data()), indexBuffer.data(), deviceIndices.numBytes()));
        CUDA_ASSERT(cuMemcpyHtoD(
            CUdeviceptr(deviceWidths.data()), hostWidths.data(), deviceWidths.numBytes()));

        unsigned int flags = OPTIX_GEOMETRY_FLAG_REQUIRE_SINGLE_ANYHIT_CALL;

        OptixBuildInput input = {};
        input.type = OPTIX_BUILD_INPUT_TYPE_CURVES;
        OptixBuildInputCurveArray &curveArray = input.curveArray;
        curveArray = {};
        curveArray.curveType = OPTIX_PRIMITIVE_TYPE_ROUND_CUBIC_BSPLINE;
        curveArray.numPrimitives = totalNumSegments;
        curveArray.vertexBuffers = vertexBuffers.data();
        curveArray.numVertices = curve.GetNumVertices();
        curveArray.vertexStrideInBytes = 0;

        curveArray.widthBuffers = widthBuffers.data();
        curveArray.widthStrideInBytes = 0;
        curveArray.indexBuffer = CUdeviceptr(deviceIndices.data());
        curveArray.indexStrideInBytes = 0;
        curveArray.flag = flags;
        curveArray.primitiveIndexOffset = 0;
        curveArray.endcapFlags = OPTIX_CURVE_ENDCAP_DEFAULT;

        BuildOptixBVH(cudaDevice, deviceArena, options, input);

        hostArena.Clear();
        deviceArena.Clear();
    }

#if (OPTIX_VERSION >= 90000)
    const int templateMinDim = 2;
    const int templateMaxDim = 8;
    const int templateWidths = templateMaxDim - templateMinDim + 1;

    for (MicropolygonMesh &mesh : scene->micropolygonMeshes)
    {
        if (cudaDevice->gridClusterTemplateOutputHandles.data() == nullptr ||
            cudaDevice->gridClusterTemplateOutputHandles.size() == 0)
            continue;

        const size_t numPositions = mesh.positions.size();
        if (numPositions == 0 || mesh.grids.size() == 0)
            continue;

        MemoryView<float3> hostGridVertices = hostArena.PushArray<float3>(numPositions);
        DeviceMemoryView<float3> deviceGridVertices = deviceArena.PushArray<float3>(numPositions);

        size_t numBuildGrids = 0;
        for (size_t gridIdx = 0; gridIdx < mesh.grids.size(); gridIdx++)
        {
            const Grid &grid = mesh.grids[gridIdx];
            int w = grid.width;
            int h = grid.height;
            if (w < templateMinDim || w > templateMaxDim || h < templateMinDim ||
                h > templateMaxDim)
                continue;
            numBuildGrids++;
        }
        if (numBuildGrids == 0)
        {
            hostArena.Clear();
            deviceArena.Clear();
            continue;
        }

        MemoryView<OptixClusterAccelBuildInputTemplatesArgs> hostTemplateArgs =
            hostArena.PushArray<OptixClusterAccelBuildInputTemplatesArgs>(numBuildGrids);
        size_t buildGridIdx = 0;
        size_t totalNumVertices = 0;

        for (size_t gridIdx = 0; gridIdx < mesh.grids.size(); gridIdx++)
        {
            const Grid &grid = mesh.grids[gridIdx];
            int w = grid.width;
            int h = grid.height;
            if (w < templateMinDim || w > templateMaxDim || h < templateMinDim ||
                h > templateMaxDim)
                continue;

            const int templateIndex = (w - templateMinDim) * templateWidths + (h - templateMinDim);
            const OptixTraversableHandle templateHandle =
                cudaDevice->gridClusterTemplateOutputHandles.data()[templateIndex];

            const int numVertices = (w + 1) * (h + 1);
            const int startVertex = grid.gridIndexStart;
            YBI_ASSERT(startVertex >= 0 &&
                       (size_t)(startVertex + numVertices) <= mesh.positions.size());

            memcpy(hostGridVertices.data() + totalNumVertices,
                   mesh.positions.data() + startVertex,
                   (size_t)numVertices * sizeof(float3));

            const CUdeviceptr vertexBufferPtr =
                (CUdeviceptr)deviceGridVertices.data() + totalNumVertices * sizeof(float3);

            OptixClusterAccelBuildInputTemplatesArgs &templateArgs =
                hostTemplateArgs[buildGridIdx];
            templateArgs = {};
            templateArgs.clusterTemplate = templateHandle;
            templateArgs.clusterIdOffset = 0;
            templateArgs.vertexBuffer = vertexBufferPtr;
            buildGridIdx++;
            totalNumVertices += (size_t)numVertices;
        }

        CUDA_ASSERT(cuMemcpyHtoD(CUdeviceptr(deviceGridVertices.data()),
                                 hostGridVertices.data(),
                                 totalNumVertices * sizeof(float3)));

        DeviceMemoryView<OptixClusterAccelBuildInputTemplatesArgs> deviceTemplateArgs =
            deviceArena.PushArray<OptixClusterAccelBuildInputTemplatesArgs>(numBuildGrids);
        CUDA_ASSERT(
            cuMemcpyHtoD(CUdeviceptr(deviceTemplateArgs.data()),
                         hostTemplateArgs.data(),
                         numBuildGrids * sizeof(OptixClusterAccelBuildInputTemplatesArgs)));

        DeviceMemoryView<uint32_t> deviceArgsCount = deviceArena.PushArray<uint32_t>(1);
        const uint32_t argsCountVal = (uint32_t)numBuildGrids;
        CUDA_ASSERT(
            cuMemcpyHtoD(CUdeviceptr(deviceArgsCount.data()), &argsCountVal, sizeof(uint32_t)));

        const unsigned int maxTrianglesPerCluster = 128;
        const unsigned int maxVerticesPerCluster = 128;

        OptixClusterAccelBuildInput clusterInput = {};
        clusterInput.type = OPTIX_CLUSTER_ACCEL_BUILD_TYPE_CLUSTERS_FROM_TEMPLATES;
        clusterInput.triangles.flags = OPTIX_CLUSTER_ACCEL_BUILD_FLAG_PREFER_FAST_TRACE;
        clusterInput.triangles.maxArgCount = (unsigned int)numBuildGrids;
        clusterInput.triangles.vertexFormat = OPTIX_VERTEX_FORMAT_FLOAT3;
        clusterInput.triangles.maxTriangleCountPerArg = maxTrianglesPerCluster;
        clusterInput.triangles.maxVertexCountPerArg = maxVerticesPerCluster;

        OptixAccelBufferSizes getSizeBufferSizes = {};
        OPTIX_ASSERT(optixClusterAccelComputeMemoryUsage(cudaDevice->optixDeviceContext,
                                                         OPTIX_CLUSTER_ACCEL_BUILD_MODE_GET_SIZES,
                                                         &clusterInput,
                                                         &getSizeBufferSizes));

        OptixAccelBufferSizes explicitBufferSizes = {};
        OPTIX_ASSERT(optixClusterAccelComputeMemoryUsage(
            cudaDevice->optixDeviceContext,
            OPTIX_CLUSTER_ACCEL_BUILD_MODE_EXPLICIT_DESTINATIONS,
            &clusterInput,
            &explicitBufferSizes));

        const size_t tempSizeInBytes =
            std::max(getSizeBufferSizes.tempSizeInBytes, explicitBufferSizes.tempSizeInBytes);
        DeviceMemoryView<size_t> deviceOutputSizes = deviceArena.PushArray<size_t>(numBuildGrids);
        DeviceMemoryView<char> temp = deviceArena.PushArray<char>(tempSizeInBytes);

        OptixClusterAccelBuildModeDesc getSizeDesc = {};
        getSizeDesc.mode = OPTIX_CLUSTER_ACCEL_BUILD_MODE_GET_SIZES;
        getSizeDesc.getSize.outputSizesBuffer = CUdeviceptr(deviceOutputSizes.data());
        getSizeDesc.getSize.tempBuffer = CUdeviceptr(temp.data());
        getSizeDesc.getSize.tempBufferSizeInBytes = getSizeBufferSizes.tempSizeInBytes;

        OPTIX_ASSERT(optixClusterAccelBuild(cudaDevice->optixDeviceContext,
                                            0,
                                            &getSizeDesc,
                                            &clusterInput,
                                            CUdeviceptr(deviceTemplateArgs.data()),
                                            CUdeviceptr(deviceArgsCount.data()),
                                            sizeof(OptixClusterAccelBuildInputTemplatesArgs)));

        CUDA_ASSERT(cuStreamSynchronize(0));

        MemoryView<size_t> hostOutputSizes = hostArena.PushArray<size_t>(numBuildGrids);
        CUDA_ASSERT(cuMemcpyDtoH(hostOutputSizes.data(),
                                 CUdeviceptr(deviceOutputSizes.data()),
                                 numBuildGrids * sizeof(size_t)));

        size_t totalOutputSize = 0;
        MemoryView<size_t> hostOffsets = hostArena.PushArray<size_t>(numBuildGrids);
        for (size_t i = 0; i < numBuildGrids; i++)
        {
            hostOffsets[i] = totalOutputSize;
            totalOutputSize += hostOutputSizes[i];
        }

        YBI_ASSERT(totalOutputSize != 0);

        DeviceMemoryView<char> explicitOutput = deviceArena.PushArray<char>(totalOutputSize);
        DeviceMemoryView<size_t> deviceOutputOffsets =
            deviceArena.PushArray<size_t>(numBuildGrids);
        CUDA_ASSERT(cuMemcpyHtoD(CUdeviceptr(deviceOutputOffsets.data()),
                                 hostOffsets.data(),
                                 numBuildGrids * sizeof(size_t)));

        OptixClusterAccelBuildModeDesc explicitDesc = {};
        explicitDesc.mode = OPTIX_CLUSTER_ACCEL_BUILD_MODE_EXPLICIT_DESTINATIONS;
        explicitDesc.explicitDest.tempBuffer = CUdeviceptr(temp.data());
        explicitDesc.explicitDest.tempBufferSizeInBytes = explicitBufferSizes.tempSizeInBytes;
        explicitDesc.explicitDest.outputBuffer = CUdeviceptr(explicitOutput.data());
        explicitDesc.explicitDest.outputOffsetsBuffer = CUdeviceptr(deviceOutputOffsets.data());
        explicitDesc.explicitDest.outputSizesBuffer = CUdeviceptr(deviceOutputSizes.data());

        OPTIX_ASSERT(optixClusterAccelBuild(cudaDevice->optixDeviceContext,
                                            0,
                                            &explicitDesc,
                                            &clusterInput,
                                            CUdeviceptr(deviceTemplateArgs.data()),
                                            CUdeviceptr(deviceArgsCount.data()),
                                            sizeof(OptixClusterAccelBuildInputTemplatesArgs)));

        CUDA_ASSERT(cuStreamSynchronize(0));

        hostArena.Clear();
        deviceArena.Clear();
    }
#else
    YBI_ERROR(scene->micropolygonMeshes.size() == 0,
              "Direct rendering of micropolygon meshes are not supported on this hardware.");
#endif

#if 0
    Array<OptixBuildInput> optixBuildInputs;


    for (int collectionIndex = 0; collectionIndex < scene->primitiveCollections.size();
         collectionIndex++)
    {
        int numPrimitives = scene->GetNumPrimitives(collectionIndex);
        optixBuildInputs.Resize(numPrimitives);
        int index = 0;

        // TODO: either enforce collections having the same number of keys,
        // or split collections here if they don't
        uint32_t numMotionKeys = 1;

        OptixAccelBuildOptions options = {};
        options.buildFlags =
            OPTIX_BUILD_FLAG_ALLOW_COMPACTION | OPTIX_BUILD_FLAG_PREFER_FAST_TRACE;
        options.operation = OPTIX_BUILD_OPERATION_BUILD;
        options.motionOptions.numKeys = numMotionKeys;
        options.motionOptions.flags = 0;
        options.motionOptions.timeBegin = 0.f;
        options.motionOptions.timeEnd = 1.f;

        for (const Primitive &primitive : scene->GetPrimitivesInCollection(collectionIndex))
        {
            switch (primitive.primitiveType)
            {
                case PRIMITIVE_TYPE_TRIANGLES:
                {
                    optixBuildInputs[index] = GetOptiXTriangleBuildInput(
                        cudaDevice, memoryArena, scene->meshes[primitive.index], numMotionKeys);
                }
                break;
                case PRIMITIVE_TYPE_CURVES:
                {
                }
                break;
                case PRIMITIVE_TYPE_INSTANCES:
                {
                }
                break;
            }
            OptixBuildInput buildInput = {};
            // buildInputs
        }
    }
#endif

#if 0
#if (OPTIX_VERSION >= 90000)
    if (bvh->flags & BVHFlags::USE_CLUSTERS)
    {
        ClusterAccelerationStructureLimits limits;
        OPTIX_ASSERT(
            optixDeviceContextGetProperty(cudaDevice->optixDeviceContext,
                                          OPTIX_DEVICE_PROPERTY_LIMIT_MAX_CLUSTER_TRIANGLES,
                                          &limits.maxTrianglesPerCluster,
                                          sizeof(unsigned int)));
        OPTIX_ASSERT(
            optixDeviceContextGetProperty(cudaDevice->optixDeviceContext,
                                          OPTIX_DEVICE_PROPERTY_LIMIT_MAX_CLUSTER_VERTICES,
                                          &limits.maxVerticesPerCluster,
                                          sizeof(unsigned int)));
        printf("limits: %u %u\n", limits.maxVerticesPerCluster, limits.maxTrianglesPerCluster);

        OptixClusterAccelBuildInputTrianglesArgs *args;
        uint32_t *argsCount;
        uint32_t *outputSizes;
        void *tempBuffer;
        size_t tempBufferSize           = 0;
        uint32_t numClusters            = 0;
        uint32_t maxTrianglesPerCluster = 128;
        uint32_t maxVerticesPerCluster  = 256;

        OptixClusterAccelBuildModeDesc buildModeDesc = {};
        buildModeDesc.mode                          = OPTIX_CLUSTER_ACCEL_BUILD_MODE_GET_SIZES;
        buildModeDesc.getSize.outputSizesBuffer     = CUdeviceptr(outputSizes);
        buildModeDesc.getSize.tempBuffer            = CUdeviceptr(tempBuffer);
        buildModeDesc.getSize.tempBufferSizeInBytes = tempBufferSize;

        OptixClusterAccelBuildInput buildInput = {};
        buildInput.type            = OPTIX_CLUSTER_ACCEL_BUILD_TYPE_CLUSTERS_FROM_TRIANGLES;
        buildInput.triangles.flags = OPTIX_CLUSTER_ACCEL_BUILD_FLAG_PREFER_FAST_TRACE;
        buildInput.triangles.maxArgCount            = numClusters;
        buildInput.triangles.vertexFormat           = OPTIX_VERTEX_FORMAT_FLOAT3;
        buildInput.triangles.maxTriangleCountPerArg = maxTrianglesPerCluster;
        buildInput.triangles.maxVertexCountPerArg   = maxVerticesPerCluster;

        OPTIX_ASSERT(optixClusterAccelBuild(cudaDevice->optixDeviceContext,
                                            0,
                                            &buildModeDesc,
                                            &buildInput,
                                            CUdeviceptr(args),
                                            CUdeviceptr(argsCount),
                                            0));
    }
#else
    YBI_ERROR(bvh->flags & BVHFlags::USE_CLUSTERS, "Cluster Acceleration Structures are not supported. Please update OptiX.");
#endif
#endif

    CUDA_ASSERT(cuCtxPopCurrent(0));
}

#endif
#endif

YBI_NAMESPACE_END
