#include "device/device.h"

#include <cstdio>

#if defined(WITH_CUDA) && defined(WITH_OPTIX)
#include "device/cuda_device.h"
#endif

YBI_NAMESPACE_BEGIN

Device *Device::CreateDevice(DeviceKind kind, std::unique_ptr<Device> &storage)
{
    if (kind == DeviceKind::CPU)
    {
        fprintf(stderr, "CPU device backend not implemented yet.\n");
        return nullptr;
    }

#if defined(WITH_CUDA) && defined(WITH_OPTIX)
    storage = std::make_unique<CUDADevice>();
    return storage.get();
#else
    (void)storage;
    fprintf(stderr, "GPU backend requires WITH_CUDA and WITH_OPTIX.\n");
    return nullptr;
#endif
}

YBI_NAMESPACE_END

