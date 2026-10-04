#include "render/wavefront_integrator.h"
#include "device/device.h"
#include "kernel/launch/wavefront_params.h"

namespace ybi
{

WavefrontIntegrator::WavefrontIntegrator(Device *device) : device(device)
{
    if (!device)
        return;

#if defined(YBI_COUNT_KERNELS_KERNEL_PATH)
    DeviceKernelDesc desc;
    desc.filename = YBI_COUNT_KERNELS_KERNEL_PATH;
    desc.entryName = "CountKernels";
    countKernels = device->CreateDeviceKernel(desc);
    if (countKernels)
        kernelCountsBuffer = device->Alloc(sizeof(KernelCounts));
#endif

#if defined(YBI_SELECT_SAMPLE_INDICES_KERNEL_PATH)
    DeviceKernelDesc selectDesc;
    selectDesc.filename = YBI_SELECT_SAMPLE_INDICES_KERNEL_PATH;
    selectDesc.entryName = "SelectSampleIndices";
    selectSampleIndices = device->CreateDeviceKernel(selectDesc);
    if (selectSampleIndices)
        numSampleIndicesBuffer = device->Alloc(sizeof(uint32_t));
#endif
}

WavefrontIntegrator::~WavefrontIntegrator()
{
    if (kernelCountsBuffer)
        device->Free(kernelCountsBuffer);
    if (numSampleIndicesBuffer)
        device->Free(numSampleIndicesBuffer);
}

bool WavefrontIntegrator::RenderSamples(const WavefrontParams &params)
{
    // uint32_t numSampleStates = 0;
    // device->CopyToHost(&numSampleStates, params.numSampleStates, sizeof(numSampleStates));
    // if (numSampleStates == 0)
    //     return true;
    // if (numSampleStates > params.maxNumSampleState || !params.sampleState.kernelTypes)
    //     return false;

    const KernelType *kernelTypes = params.sampleState.kernelTypes;
    KernelCounts *deviceCounts = kernelCountsBuffer.As<KernelCounts>().Get();
    void *arguments[] = {&kernelTypes, &numSampleStates, &deviceCounts};

    DeviceKernelLaunchDesc launch;
    launch.blockX = 256;
    launch.gridX = (numSampleStates - 1) / launch.blockX + 1;
    launch.arguments = arguments;

    while (true)
    {
        device->Memzero(kernelCountsBuffer.Data(), sizeof(KernelCounts));
        device->ExecuteKernel(*countKernels, launch);

        KernelCounts counts;
        device->CopyToHost(&counts, kernelCountsBuffer.Data(), sizeof(counts));
        if (counts.counts[KERNEL_FINISHED] == numSampleStates)
            return true;

        if (!params.sampleIndices)
            return false;

        KernelType nextKernel = KERNEL_FINISHED;
        for (int type = KERNEL_INTERSECT; type < KERNEL_FINISHED; ++type)
        {
            if (counts.counts[type] > 0)
            {
                nextKernel = static_cast<KernelType>(type);
                break;
            }
        }
        if (nextKernel == KERNEL_FINISHED)
            return false;

        WavefrontParams dispatchParams = params;
        void *selectArguments[] = {&kernelTypes,
                                   &numSampleStates,
                                   &nextKernel,
                                   &dispatchParams.sampleIndices,
                                   numSampleIndicesBuffer.KernelArg()};
        DeviceKernelLaunchDesc selectLaunch = launch;
        selectLaunch.arguments = selectArguments;

        device->Memzero(numSampleIndicesBuffer.Data(), sizeof(uint32_t));
        device->ExecuteKernel(*selectSampleIndices, selectLaunch);
        device->CopyToHost(&dispatchParams.numSampleIndices,
                           numSampleIndicesBuffer.Data(),
                           sizeof(dispatchParams.numSampleIndices));

        // TODO: upload dispatchParams and dispatch nextKernel before recounting.
        return false;
    }
}

} // namespace ybi
