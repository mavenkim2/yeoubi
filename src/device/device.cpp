#include "device/device.h" 
#include "device/cpu_device.h"
#if defined(WITH_CUDA) && defined(WITH_OPTIX)
#include "device/cuda_device.h"
#endif

namespace ybi 
{

std::unique_ptr<Device> Device::CreateDevice(DeviceKind kind, bool createRendererPipeline)
{
    if (kind == DeviceKind::CPU)
    {
        return std::make_unique<CPUDevice>(createRendererPipeline);
    }

    if (kind == DeviceKind::GPU)
    {
#if defined(WITH_CUDA) && defined(WITH_OPTIX)
        return std::make_unique<CUDADevice>(createRendererPipeline);
#else
        fprintf(stderr, "GPU backend requires WITH_CUDA and WITH_OPTIX.\n");
        return nullptr;
#endif
    }

    fprintf(stderr, "Unknown device kind.\n");
    return nullptr;
}

}
