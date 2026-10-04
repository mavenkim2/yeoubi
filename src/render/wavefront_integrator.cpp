#include "render/wavefront_integrator.h"
#include "core/time.h"
#include "device/device.h"
#include "kernel/launch/wavefront_params.h"
#include <cstdio>
#include <map>
#include <nvtx3/nvtx3.hpp>
#include <string>
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
    kernelCountsBuffer = device->Alloc(sizeof(KernelCounts));

    DeviceKernelDesc selectDesc;
    selectDesc.filename = YBI_SELECT_SAMPLE_INDICES_KERNEL_PATH;
    selectDesc.entryName = "SelectSampleIndices";
    selectSampleIndices = device->CreateDeviceKernel(selectDesc);
    numSampleIndicesBuffer = device->Alloc(sizeof(uint32_t));

    DeviceKernelDesc compactionDesc;
    compactionDesc.filename = YBI_STREAM_COMPACTION_KERNEL_PATH;
    compactionDesc.entryName = "StreamCompaction";
    streamCompaction = device->CreateDeviceKernel(compactionDesc);

    AllocateSampleState(sampleState);
    numSampleStatesBuffer = device->Alloc(sizeof(uint32_t));
    sampleIndicesBuffer = device->Alloc(sizeof(uint32_t) * batchSize);
    kernelParamsBuffer = device->Alloc(sizeof(WavefrontParams));
}

WavefrontIntegrator::~WavefrontIntegrator()
{
    FreeSampleState();
    if (kernelCountsBuffer)
        device->Free(kernelCountsBuffer);
    if (numSampleIndicesBuffer)
        device->Free(numSampleIndicesBuffer);
    if (numSampleStatesBuffer)
        device->Free(numSampleStatesBuffer);
    if (sampleIndicesBuffer)
        device->Free(sampleIndicesBuffer);
    if (kernelParamsBuffer)
        device->Free(kernelParamsBuffer);
}

void WavefrontIntegrator::AllocateSampleState(WavefrontSampleState &state)
{
    const auto allocate = [this](auto *&ptr, size_t elementsPerSample = 1) {
        using T = std::remove_pointer_t<std::remove_reference_t<decltype(ptr)>>;
        DeviceBuffer buffer = device->Alloc(sizeof(T) * batchSize * elementsPerSample);
        sampleStateBuffers.push_back(buffer);
        ptr = buffer.As<T>().Get();
    };

    allocate(state.pixelIndex);
    allocate(state.rng);
    allocate(state.rayCone);
    allocate(state.radiance);
    allocate(state.throughput);
    allocate(state.previousSurfacePosition);
    allocate(state.previousBsdfPdf);
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
    for (DeviceBuffer &buffer : sampleStateBuffers)
        device->Free(buffer);
    sampleStateBuffers.clear();
    sampleState = {};
    compactionState = {};
}

bool WavefrontIntegrator::RenderSamples(const KernelParams &kernelParams)
{
    const nvtx3::scoped_range renderRange{"RenderSamples"};
    std::map<std::string, double> timings;
    const auto finish = [&timings](bool success) {
        for (const auto &[name, milliseconds] : timings)
            std::printf("%s: %.3f ms total\n", name.c_str(), milliseconds);
        return success;
    };
    const auto dispatchKernel = [this, &timings](KernelType type,
                                                     uint32_t numWorkItems,
                                                     uint32_t blockSize) {
        static constexpr const char *names[] = {
            "GenerateCameraRays", "Intersect", "ShadeHit", "DirectLighting", "IndirectLighting"};
        const nvtx3::scoped_range range{names[type]};
        float milliseconds = 0.0f;
        const bool success = device->DispatchKernel(
            type, numWorkItems, blockSize, kernelParamsBuffer, &milliseconds);
        timings[names[type]] += milliseconds;
        return success;
    };
    const auto executeKernel = [this, &timings](const char *name,
                                                    DeviceKernel &kernel,
                                                    const DeviceKernelLaunchDesc &launch) {
        const nvtx3::scoped_range range{name};
        float milliseconds = 0.0f;
        device->ExecuteKernel(kernel, launch, &milliseconds);
        timings[name] += milliseconds;
    };
    const auto memzero = [this, &timings](void *dst, size_t size) {
        const nvtx3::scoped_range range{"Memzero"};
        const Timer timer = StartTimer();
        device->Memzero(dst, size);
        const double milliseconds = StopTimer(timer);
        timings["Memzero (wall)"] += milliseconds;
    };
    const auto copyToHost = [this, &timings](void *dst, void *src, size_t size) {
        const nvtx3::scoped_range range{"CopyToHost"};
        const Timer timer = StartTimer();
        device->CopyToHost(dst, src, size);
        const double milliseconds = StopTimer(timer);
        timings["CopyToHost (wall)"] += milliseconds;
    };
    const auto copyToDevice = [this, &timings](void *dst, void *src, size_t size) {
        const nvtx3::scoped_range range{"CopyToDevice"};
        const Timer timer = StartTimer();
        device->CopyToDevice(dst, src, size);
        const double milliseconds = StopTimer(timer);
        timings["CopyToDevice (wall)"] += milliseconds;
    };

    WavefrontParams params{};
    static_cast<KernelParams &>(params) = kernelParams;
    params.sampleState = sampleState;
    params.numSampleStates = numSampleStatesBuffer.As<uint32_t>().Get();
    params.maxNumSampleState = batchSize;
    params.sampleIndices = sampleIndicesBuffer.As<uint32_t>().Get();

    const uint64_t numPixels = uint64_t(params.imageWidth) * params.imageHeight;
    uint32_t numSampleStates = static_cast<uint32_t>(Min<uint64_t>(batchSize, numPixels));

    const KernelType *kernelTypes = params.sampleState.kernelTypes;
    KernelCounts *deviceCounts = kernelCountsBuffer.As<KernelCounts>().Get();
    void *arguments[] = {&kernelTypes, &numSampleStates, &deviceCounts};

    DeviceKernelLaunchDesc launch;
    launch.blockX = 256;
    launch.gridX = (numSampleStates - 1) / launch.blockX + 1;
    launch.arguments = arguments;

    WavefrontParams initialParams = params;
    initialParams.sampleIndices = nullptr;
    initialParams.numSampleIndices = numSampleStates;

    uint64_t nextWorkIndex = 0;
    initialParams.cameraWorkOffset = nextWorkIndex;
    initialParams.cameraSampleStateOffset = 0;
    copyToDevice(params.numSampleStates, &numSampleStates, sizeof(numSampleStates));
    copyToDevice(kernelParamsBuffer.Data(), &initialParams, sizeof(initialParams));
    if (!dispatchKernel(KERNEL_GENERATE_CAMERA_RAYS, numSampleStates, launch.blockX))
        return finish(false);
    nextWorkIndex += numSampleStates;

    if (!dispatchKernel(KERNEL_INTERSECT, numSampleStates, launch.blockX))
        return finish(false);

    while (true)
    {
        memzero(kernelCountsBuffer.Data(), sizeof(KernelCounts));
        executeKernel("CountKernels", *countKernels, launch);

        KernelCounts counts;
        copyToHost(&counts, kernelCountsBuffer.Data(), sizeof(counts));
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
            if (!compactionState.kernelTypes)
                AllocateSampleState(compactionState);

            void *compactionArguments[] = {&sampleState,
                                           &compactionState,
                                           &numSampleStates,
                                           numSampleStatesBuffer.KernelArg()};
            DeviceKernelLaunchDesc compactionLaunch = launch;
            compactionLaunch.arguments = compactionArguments;
            memzero(numSampleStatesBuffer.Data(), sizeof(uint32_t));
            executeKernel("StreamCompaction", *streamCompaction, compactionLaunch);
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
            if (!dispatchKernel(KERNEL_GENERATE_CAMERA_RAYS, numCameraSamples, launch.blockX))
                return finish(false);

            nextWorkIndex += numCameraSamples;
            numSampleStates += numCameraSamples;
            counts.counts[KERNEL_INTERSECT] += numCameraSamples;
            copyToDevice(
                params.numSampleStates, &numSampleStates, sizeof(numSampleStates));
            launch.gridX = (numSampleStates - 1) / launch.blockX + 1;
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

        memzero(numSampleIndicesBuffer.Data(), sizeof(uint32_t));
        executeKernel("SelectSampleIndices", *selectSampleIndices, selectLaunch);

        copyToDevice(kernelParamsBuffer.Data(), &dispatchParams, sizeof(dispatchParams));
        if (!dispatchKernel(nextKernel, dispatchParams.numSampleIndices, launch.blockX))
            return finish(false);
    }
}

} // namespace ybi
