#include "kernel/integrator/generate_camera_ray.h"

namespace ybi
{

extern "C" __global__ void GenerateCameraRays(WavefrontParams params)
{
    const uint32_t x = blockDim.x * blockIdx.x + threadIdx.x;
    const uint32_t y = blockDim.y * blockIdx.y + threadIdx.y;
    const Vec2u launchDimensions(static_cast<uint32_t>(params.imageWidth),
                                 static_cast<uint32_t>(params.imageHeight));
    if (x >= launchDimensions.x || y >= launchDimensions.y)
        return;

    GenerateCameraRay(params, Vec2u(x, y), launchDimensions);
}

} // namespace ybi
