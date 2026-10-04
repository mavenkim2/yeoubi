#pragma once

#include "core/assert.h"
#include "core/vec2u.h"
#include "kernel/launch/wavefront_params.h"
#include "kernel/sampling/rng.h"

namespace ybi
{

// Generate one sample per launch pixel. The launch must complete before tracing
// consumes the sample count or state buffers.
YBI_GPU void GenerateCameraRay(const WavefrontParams &params,
                               const Vec2u &launchIndex,
                               const Vec2u &launchDimensions)
{

    const uint64_t sampleCount = uint64_t(launchDimensions.x) * launchDimensions.y;
    YBI_ASSERT(sampleCount <= params.maxNumSampleState);

    const uint32_t sampleIndex = launchIndex.x + launchIndex.y * launchDimensions.x;
    const Vec2u raster(launchIndex.x + params.pixelOffsetX, launchIndex.y + params.pixelOffsetY);
    RNG rng = RNG::Init(raster, params.sampleOffset);
    const float pixelX = float(raster.x) + rng.Uniform();
    const float pixelY = float(raster.y) + rng.Uniform();
    Vec4f cameraPos = params.cameraFromRaster * Vec4f(Vec3f(pixelX, pixelY, 0.f), 1.f);
    cameraPos.xyz /= cameraPos.w;

    const Vec3f rayOrigin = params.renderFromCamera.TransformPoint(Vec3f(0.f));
    const Vec3f rayDirection = Normalize(params.renderFromCamera.TransformVector(cameraPos.xyz));
    RayCone rayCone{};
    rayCone.spreadAngle = atanf(2.f * tanf(0.5f * params.fovRadians) / float(params.imageHeight));

    const WavefrontSampleState &state = params.sampleState;
    state.WritePixelIndex(sampleIndex, sampleIndex);
    state.WriteRNG(sampleIndex, rng);
    state.WriteRayCone(sampleIndex, rayCone);
    state.WriteRadiance(sampleIndex, Vec3f(0.f));
    state.WriteThroughput(sampleIndex, Vec3f(1.f));
    state.WritePreviousSurfacePosition(sampleIndex, Vec3f(0.f));
    state.WritePreviousBsdfPdf(sampleIndex, 0.f);
    state.WriteRayOrigin(sampleIndex, rayOrigin);
    state.WriteRayDirection(sampleIndex, rayDirection);
    state.WriteDepth(sampleIndex, 0u);
    state.WriteKernelType(sampleIndex, KERNEL_INTERSECT);
}

} // namespace ybi
