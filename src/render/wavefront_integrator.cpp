#include "render/wavefront_integrator.h"
#include "core/math.h"
#include "device/device.h"
#include "kernel/guiding/guiding.h"
#include "kernel/guiding/kd_tree.h"
#include "kernel/launch/wavefront_params.h"
#include "render/path_guiding/path_guiding.h"
#include <cstring>
#include <type_traits>
#include <utility>

namespace ybi
{

WavefrontIntegrator::WavefrontIntegrator(Device *device) : device(device)
{
    if (!device)
        return;

    DeviceKernelDesc desc;
    desc.filename = YBI_COUNT_KERNELS_KERNEL_PATH;
    desc.entryName = "CountKernels";
    countKernels = device->CreateDeviceKernel(desc);
    kernelCountsBuffer = device->Alloc<KernelCounts>(1);
    hostKernelCountsBuffer = device->AllocHost<KernelCounts>(1);
    // A refill can upload camera params, the new state count, and dispatch params
    // before the next readback. Each upload keeps its own bytes until that sync.
    hostUploadBuffer = device->AllocHost(3 * sizeof(WavefrontParams));

    DeviceKernelDesc selectDesc;
    selectDesc.filename = YBI_SELECT_SAMPLE_INDICES_KERNEL_PATH;
    selectDesc.entryName = "SelectSampleIndices";
    selectSampleIndices = device->CreateDeviceKernel(selectDesc);
    numSampleIndicesBuffer = device->Alloc<uint32_t>(1);

    DeviceKernelDesc compactionDesc;
    compactionDesc.filename = YBI_STREAM_COMPACTION_KERNEL_PATH;
    compactionDesc.entryName = "StreamCompaction";
    streamCompaction = device->CreateDeviceKernel(compactionDesc);

    AllocateSampleState(sampleState);
    AllocateSampleState(compactionState);
    numSampleStatesBuffer = device->Alloc<uint32_t>(1);
    sampleIndicesBuffer = device->Alloc<uint32_t>(batchSize);
    kernelParamsBuffer = device->Alloc<WavefrontParams>(1);
#if defined(YBI_USE_PATH_GUIDING)
    radianceRecordAllocationCursorBuffer = device->Alloc<uint32_t>(1);
    radianceRecordRetiredCountBuffer = device->Alloc<uint32_t>(1);
    pathGuidingTrainingSampleCountBuffer = device->Alloc<uint32_t>(1);
#endif
}

WavefrontIntegrator::~WavefrontIntegrator()
{
    FreeSampleState();
    if (kernelCountsBuffer)
        device->Free(kernelCountsBuffer);
    if (hostKernelCountsBuffer)
        device->FreeHost(hostKernelCountsBuffer);
    if (hostUploadBuffer)
        device->FreeHost(hostUploadBuffer);
    if (numSampleIndicesBuffer)
        device->Free(numSampleIndicesBuffer);
    if (numSampleStatesBuffer)
        device->Free(numSampleStatesBuffer);
    if (sampleIndicesBuffer)
        device->Free(sampleIndicesBuffer);
    if (kernelParamsBuffer)
        device->Free(kernelParamsBuffer);
#if defined(YBI_USE_PATH_GUIDING)
    if (radianceRecordAllocationCursorBuffer)
        device->Free(radianceRecordAllocationCursorBuffer);
    if (radianceRecordRetiredCountBuffer)
        device->Free(radianceRecordRetiredCountBuffer);
    if (pathGuidingTrainingSamplesBuffer)
        device->Free(pathGuidingTrainingSamplesBuffer);
    if (pathGuidingTrainingSampleCountBuffer)
        device->Free(pathGuidingTrainingSampleCountBuffer);
#endif
}

void WavefrontIntegrator::AllocateSampleState(WavefrontSampleState &state)
{
    const auto allocate = [this](auto *&ptr, size_t elementsPerSample = 1) {
        using T = std::remove_pointer_t<std::remove_reference_t<decltype(ptr)>>;
        DeviceBuffer<T> buffer = device->Alloc<T>(batchSize * elementsPerSample);
        sampleStateBuffers.push_back(buffer.template As<std::byte>());
        ptr = buffer.Data();
    };

    allocate(state.pixelIndex);
    allocate(state.rng);
    allocate(state.rayCone);
    allocate(state.radiance);
#if defined(YBI_USE_PATH_GUIDING)
    allocate(state.radianceRecordIndex);
#endif
    allocate(state.throughput);
    allocate(state.previousSurfacePosition);
    allocate(state.previousBsdfPdf);
    allocate(state.previousBsdfFlags);
    allocate(state.etaScale);
    allocate(state.rayOrigin);
    allocate(state.rayDirection);
    allocate(state.depth);
    allocate(state.kernelTypes);

    allocate(state.t);
    allocate(state.bary);
    allocate(state.objectToWorld);
    allocate(state.primIndex);
    allocate(state.clusterId);
    allocate(state.objectId);
    allocate(state.type);

    allocate(state.materialIndex);
    allocate(state.position);
    allocate(state.geometricNormal);
    allocate(state.shadingNormal);
    allocate(state.dpdu);
    allocate(state.uv);
    allocate(state.triangleArea);
    allocate(state.uvArea);
    allocate(state.closureData, MAX_CLOSURE_SIZE);
}

void WavefrontIntegrator::FreeSampleState()
{
    for (DeviceBuffer<std::byte> &buffer : sampleStateBuffers)
        device->Free(buffer);
    sampleStateBuffers.clear();
    sampleState = {};
    compactionState = {};
}

bool WavefrontIntegrator::RenderSamples(const KernelParams &kernelParams)
{
    size_t uploadOffset = 0;
    const auto copyToDevice = [this, &uploadOffset](void *dst, void *src, size_t size) {
        YBI_ASSERT(uploadOffset + size <= hostUploadBuffer.SizeInBytes());
        auto *staging = hostUploadBuffer.Data() + uploadOffset;
        std::memcpy(staging, src, size);
        device->CopyToDeviceAsync(dst, staging, size);
        uploadOffset += size;
    };
    const auto finish = [this](bool success) {
        // Failed dispatches can leave uploads in flight. Drain them before the
        // next frame reuses the pinned staging buffer.
        if (!success)
            device->SyncDeviceStream();
        return success;
    };

    WavefrontParams params{};
    static_cast<KernelParams &>(params) = kernelParams;
#if defined(YBI_USE_PATH_GUIDING)
    if (pathGuidingTrainingSamplesBuffer.Count() < params.radianceRecorder.availableCount)
    {
        if (pathGuidingTrainingSamplesBuffer)
            device->Free(pathGuidingTrainingSamplesBuffer);
        pathGuidingTrainingSamplesBuffer =
            device->Alloc<PathGuidingTrainingSample>(params.radianceRecorder.availableCount);
    }
    params.radianceRecorder.allocationCursor = radianceRecordAllocationCursorBuffer.Data();
    params.radianceRecorder.retiredListCount = radianceRecordRetiredCountBuffer.Data();
    params.radianceRecorder.trainingSamples = pathGuidingTrainingSamplesBuffer.Data();
    params.radianceRecorder.trainingSampleCount = pathGuidingTrainingSampleCountBuffer.Data();
    params.radianceRecorder.trainingSamplesMax =
        static_cast<uint32_t>(pathGuidingTrainingSamplesBuffer.Count());
    device->MemzeroAsync(radianceRecordAllocationCursorBuffer.Data(), sizeof(uint32_t));
    device->MemzeroAsync(radianceRecordRetiredCountBuffer.Data(), sizeof(uint32_t));
    device->MemzeroAsync(pathGuidingTrainingSampleCountBuffer.Data(), sizeof(uint32_t));
#endif
    params.sampleState = sampleState;
    params.numSampleStates = numSampleStatesBuffer.Data();
    params.maxNumSampleState = batchSize;
    params.sampleIndices = sampleIndicesBuffer.Data();

    const uint64_t numPixels = uint64_t(params.imageWidth) * params.imageHeight;
    uint32_t numSampleStates = static_cast<uint32_t>(Min<uint64_t>(batchSize, numPixels));

    const KernelType *kernelTypes = params.sampleState.kernelTypes;
    KernelCounts *deviceCounts = kernelCountsBuffer.Data();
    void *arguments[] = {&kernelTypes, &numSampleStates, &deviceCounts};

    DeviceKernelLaunchDesc launch;
    launch.blockX = 256;
    launch.gridX = CeilDiv(numSampleStates, launch.blockX);
    launch.arguments = arguments;

    WavefrontParams initialParams = params;
    initialParams.sampleIndices = nullptr;
    initialParams.numSampleIndices = numSampleStates;

    uint64_t nextWorkIndex = 0;
    initialParams.cameraWorkOffset = nextWorkIndex;
    initialParams.cameraSampleStateOffset = 0;
    copyToDevice(params.numSampleStates, &numSampleStates, sizeof(numSampleStates));
    copyToDevice(kernelParamsBuffer.Data(), &initialParams, sizeof(initialParams));
    if (!device->DispatchKernel(
            KERNEL_GENERATE_CAMERA_RAYS, numSampleStates, launch.blockX, kernelParamsBuffer))
        return finish(false);
    nextWorkIndex += numSampleStates;

    if (!device->DispatchKernel(
            KERNEL_INTERSECT, numSampleStates, launch.blockX, kernelParamsBuffer))
        return finish(false);

    while (true)
    {
        device->MemzeroAsync(kernelCountsBuffer.Data(), sizeof(KernelCounts));
        device->ExecuteKernel(*countKernels, launch);

        KernelCounts &counts = *hostKernelCountsBuffer.Data();
        device->CopyToHostAsync(&counts, kernelCountsBuffer.Data(), sizeof(counts));
        device->SyncDeviceStream();
        uploadOffset = 0;
        const uint32_t numRemainingSamples = numSampleStates - counts.counts[KERNEL_FINISHED];
        if (numRemainingSamples == 0 && nextWorkIndex == numPixels)
            return finish(true);

        KernelType nextKernel = KERNEL_INTERSECT;
        uint32_t mostWork = 0;
        for (int type = KERNEL_INTERSECT; type < KERNEL_FINISHED; ++type)
        {
            if (counts.counts[type] > mostWork)
            {
                nextKernel = static_cast<KernelType>(type);
                mostWork = counts.counts[type];
            }
        }

        if (nextKernel == KERNEL_INTERSECT && numRemainingSamples < batchSize / 4 &&
            nextWorkIndex < numPixels)
        {
            void *compactionArguments[] = {&sampleState,
                                           &compactionState,
                                           &numSampleStates,
                                           numSampleStatesBuffer.KernelArg()};
            DeviceKernelLaunchDesc compactionLaunch = launch;
            compactionLaunch.arguments = compactionArguments;
            device->MemzeroAsync(numSampleStatesBuffer.Data(), sizeof(uint32_t));
            device->ExecuteKernel(*streamCompaction, compactionLaunch);
            numSampleStates = numRemainingSamples;

            std::swap(sampleState, compactionState);
            params.sampleState = sampleState;
            kernelTypes = sampleState.kernelTypes;

            const uint32_t numCameraSamples = static_cast<uint32_t>(
                Min<uint64_t>(batchSize - numSampleStates, numPixels - nextWorkIndex));
            WavefrontParams cameraParams = params;
            cameraParams.cameraWorkOffset = nextWorkIndex;
            cameraParams.cameraSampleStateOffset = numSampleStates;
            cameraParams.numSampleIndices = numCameraSamples;
            copyToDevice(kernelParamsBuffer.Data(), &cameraParams, sizeof(cameraParams));
            if (!device->DispatchKernel(KERNEL_GENERATE_CAMERA_RAYS,
                                        numCameraSamples,
                                        launch.blockX,
                                        kernelParamsBuffer))
                return finish(false);

            nextWorkIndex += numCameraSamples;
            numSampleStates += numCameraSamples;
            counts.counts[KERNEL_INTERSECT] += numCameraSamples;
            copyToDevice(params.numSampleStates, &numSampleStates, sizeof(numSampleStates));
            launch.gridX = CeilDiv(numSampleStates, launch.blockX);
        }

        WavefrontParams dispatchParams = params;
        dispatchParams.numSampleIndices = counts.counts[nextKernel];
        void *selectArguments[] = {&kernelTypes,
                                   &numSampleStates,
                                   &nextKernel,
                                   &dispatchParams.sampleIndices,
                                   numSampleIndicesBuffer.KernelArg()};
        DeviceKernelLaunchDesc selectLaunch = launch;
        selectLaunch.arguments = selectArguments;

        device->MemzeroAsync(numSampleIndicesBuffer.Data(), sizeof(uint32_t));
        device->ExecuteKernel(*selectSampleIndices, selectLaunch);

        copyToDevice(kernelParamsBuffer.Data(), &dispatchParams, sizeof(dispatchParams));
        if (!device->DispatchKernel(
                nextKernel, dispatchParams.numSampleIndices, launch.blockX, kernelParamsBuffer))
            return finish(false);
    }
}

} // namespace ybi
