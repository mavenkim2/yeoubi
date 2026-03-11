#include "device/device.h"

#include <cstdio>

#if defined(WITH_CUDA) && defined(WITH_OPTIX)
#include "device/cuda_device.h"
#endif
#if defined(WITH_EMBREE)
#include "device/cpu_device.h"
#endif

namespace ybi
{

Device *Device::CreateDevice(DeviceKind kind, std::unique_ptr<Device> &storage)
{
    if (kind == DeviceKind::CPU)
    {
#if defined(WITH_EMBREE)
        storage = std::make_unique<CPUDevice>();
        return storage.get();
#else
        (void)storage;
        fprintf(stderr, "CPU backend requires WITH_EMBREE.\n");
        return nullptr;
#endif
    }

    if (kind == DeviceKind::GPU)
    {
#if defined(WITH_CUDA) && defined(WITH_OPTIX)
        storage = std::make_unique<CUDADevice>();
        return storage.get();
#else
        (void)storage;
        fprintf(stderr, "GPU backend requires WITH_CUDA and WITH_OPTIX.\n");
        return nullptr;
#endif
    }

    (void)storage;
    fprintf(stderr, "Unknown device kind.\n");
    return nullptr;
}

} // namespace ybi
