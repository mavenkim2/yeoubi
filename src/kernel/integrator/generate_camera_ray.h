#pragma once

#include "core/assert.h"
#include "core/vec2u.h"
#include "kernel/launch/wavefront_params.h"
#include "kernel/sampling/rng.h"

namespace ybi
{

constexpr uint32_t CAMERA_TILE_WIDTH = 8;
constexpr uint32_t CAMERA_TILE_HEIGHT = 8;

// Work indices cover exactly imageWidth * imageHeight pixels, including clipped edge tiles.
YBI_CPU_GPU Vec2u CameraWorkIndexToPixel(uint64_t workIndex,
                                        uint32_t imageWidth,
                                        uint32_t imageHeight)
{
    YBI_ASSERT(imageWidth > 0 && imageHeight > 0);
    YBI_ASSERT(workIndex < uint64_t(imageWidth) * imageHeight);

    const uint64_t workPerTileRow = uint64_t(imageWidth) * CAMERA_TILE_HEIGHT;
    const uint32_t tileY = static_cast<uint32_t>(workIndex / workPerTileRow);
    const uint32_t tileHeight = Min(CAMERA_TILE_HEIGHT, imageHeight - tileY * CAMERA_TILE_HEIGHT);
    const uint64_t workInTileRow = workIndex % workPerTileRow;

    const uint32_t workPerFullTile = CAMERA_TILE_WIDTH * tileHeight;
    const uint32_t tileX = static_cast<uint32_t>(workInTileRow / workPerFullTile);
    const uint32_t tileWidth = Min(CAMERA_TILE_WIDTH, imageWidth - tileX * CAMERA_TILE_WIDTH);
    const uint32_t pixelInTile = static_cast<uint32_t>(workInTileRow % workPerFullTile);
    return Vec2u(tileX * CAMERA_TILE_WIDTH + pixelInTile % tileWidth,
                 tileY * CAMERA_TILE_HEIGHT + pixelInTile / tileWidth);
}

// Generate one camera sample into the supplied state slot. The launch must
// complete before tracing consumes the sample count or state buffers.
YBI_GPU void GenerateCameraRay(const WavefrontParams &params,
                               uint32_t sampleIndex,
                               uint64_t workIndex)
{
    YBI_ASSERT(sampleIndex < params.maxNumSampleState);
    YBI_ASSERT(params.imageWidth > 0 && params.imageHeight > 0);

    const WavefrontSampleState &state = params.sampleState;
    const Vec2u pixel = CameraWorkIndexToPixel(workIndex, params.imageWidth, params.imageHeight);

    const Vec2u raster(pixel.x + params.pixelOffsetX, pixel.y + params.pixelOffsetY);
    RNG rng = RNG::Init(raster, params.sampleOffset);
    const float pixelX = float(raster.x) + rng.Uniform();
    const float pixelY = float(raster.y) + rng.Uniform();
    Vec4f cameraPos = params.cameraFromRaster * Vec4f(Vec3f(pixelX, pixelY, 0.f), 1.f);
    cameraPos.xyz /= cameraPos.w;

    const Vec3f rayOrigin = params.renderFromCamera.TransformPoint(Vec3f(0.f));
    const Vec3f rayDirection = Normalize(params.renderFromCamera.TransformVector(cameraPos.xyz));
    RayCone rayCone{};
    rayCone.spreadAngle = atanf(2.f * tanf(0.5f * params.fovRadians) / float(params.imageHeight));

    state.WritePixelIndex(sampleIndex, pixel.x + pixel.y * params.imageWidth);
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
