#pragma once

#include "core/affine_transform.h"
#include <cstdint>

namespace ybi
{

// Tracks the view used by the last submitted frame, independently of swapchain images.
struct ProgressiveState
{
    AffineTransform camera = AffineTransform::Identity();
    uint64_t textureRevision = 0;
    uint32_t samples = 0;
    bool initialized = false;

    void BeginFrame(const AffineTransform &nextCamera, uint64_t nextTextureRevision,
                    uint32_t batchSize, bool resetRequested)
    {
        bool changed = !initialized || textureRevision != nextTextureRevision;
        for (int row = 0; row < 3; ++row)
            for (int column = 0; column < 4; ++column)
                changed |= camera.m[row][column] != nextCamera.m[row][column];
        // Keep sample counts exactly representable as floats in the running mean.
        if (changed || resetRequested || samples > (1u << 24) - batchSize)
            samples = 0;
        camera = nextCamera;
        textureRevision = nextTextureRevision;
        initialized = true;
    }
};

} // namespace ybi
