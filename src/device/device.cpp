#include "device/device.h" 
#include "device/cpu_device.h"
#include <chrono>
#if defined(WITH_CUDA) && defined(WITH_OPTIX)
#include "device/cuda_device.h"
#endif

namespace ybi 
{

namespace
{
struct HostEventTimer : DeviceEventTimer
{
    std::chrono::steady_clock::time_point start = std::chrono::steady_clock::now();
    std::chrono::steady_clock::time_point stop;
    bool stopped = false;

    float ElapsedMilliseconds() const override
    {
        YBI_ASSERT(stopped);
        return std::chrono::duration<float, std::milli>(stop - start).count();
    }
};
}

std::unique_ptr<DeviceEventTimer> Device::StartEventTimer(DeviceStream *stream)
{
    YBI_ASSERT(!stream);
    return std::make_unique<HostEventTimer>();
}

void Device::StopEventTimer(DeviceEventTimer &timer, DeviceStream *stream)
{
    YBI_ASSERT(!stream);
    auto &hostTimer = static_cast<HostEventTimer &>(timer);
    hostTimer.stop = std::chrono::steady_clock::now();
    hostTimer.stopped = true;
}

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
