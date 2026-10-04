#include "kernel/backend/optix/intersection_programs.h"
#include "kernel/integrator/primary_ray.h"
#include <optix_device.h>

namespace ybi
{

extern "C"
{
    __constant__ WavefrontParams params;
}

extern "C" __global__ void __raygen__primary()
{
    const uint3 launchIndex = optixGetLaunchIndex();
    const uint3 launchDimensions = optixGetLaunchDimensions();
    PrimaryRay(params,
               Vec2u(launchIndex.x, launchIndex.y),
               Vec2u(launchDimensions.x, launchDimensions.y));
}

} // namespace ybi
