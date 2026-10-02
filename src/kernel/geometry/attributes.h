#pragma once

#include "core/assert.h"
#include "device/cuda_compat.h"
#include "kernel/scene.h"
namespace ybi
{

YBI_CPU_GPU DeviceAttribute
GetAttribute(const DeviceScene *scene, AttributeId id, int offset, int count, int &attributeIndex)
{
    for (int index = offset; index < offset + count; index++)
    {
        const DeviceAttribute &attribute = scene->deviceAttributes[index];
        if (attribute.id == id)
        {
            attributeIndex = index;
            return attribute;
        }
    }

    YBI_LOGFATAL("Attribute with id (%u) was requested but not found\n", id);
    return {};
}

} // namespace ybi
