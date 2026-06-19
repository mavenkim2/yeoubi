#include "device/device.h" 
#include "device/cuda_device.h"

namespace ybi 
{

std::unique_ptr<Device> Device::CreateDevice(DeviceKind kind)
{
    if (kind == DeviceKind::CPU)
    {
#if defined(WITH_EMBREE)
        return std::make_unique<CPUDevice>();
#else
        fprintf(stderr, "CPU backend requires WITH_EMBREE.\n");
        return nullptr;
#endif
    }

    if (kind == DeviceKind::GPU)
    {
#if defined(WITH_CUDA) && defined(WITH_OPTIX)
        return std::make_unique<CUDADevice>();
#else
        fprintf(stderr, "GPU backend requires WITH_CUDA and WITH_OPTIX.\n");
        return nullptr;
#endif
    }

    fprintf(stderr, "Unknown device kind.\n");
    return nullptr;
}

}
