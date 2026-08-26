#pragma once

#include "core/vec2f.h"
#include "core/vec3f.h"
#include "device/device_memory.h"
#include "scene/device_cluster.h"
#include "scene/openpbr_parameters.h"
#include <cstdint>

namespace ybi
{

enum class AttributeType
{
    None,
    Vertex,
    FaceVarying,
};

struct DeviceMaterial
{
    // NOTE: = ~0u if no udims
    uint32_t udimBase;
    uint32_t udimTableOffset;

    bool valid;
};

struct DeviceMaterial2
{
#define DECLARE(name, type) int name;
    OPENPBR_PARAMETERS(DECLARE)
#undef DECLARE
};

struct Object
{
    AttributeType uvType;
    AttributeType normalType;
    uint32_t uvOffset;
    uint32_t uvIndexOffset;
    uint32_t normalOffset;
    uint32_t positionOffset;
    uint32_t triangleIndexOffset;
    uint32_t triangleMaterialOffset;
    uint32_t clusterOffset;

    int materialIndex;
};

struct DeviceScene
{
    Object *objects;
    DeviceMaterial *materials;

    Vec3f *positions;
    uint32_t *normals;
    Vec2f *uvs;
    int *uvIndices;
    int *triangleIndices;
    int *triangleMaterialIndices;

    DeviceCluster *compressedClusters;
    uint8_t *compressedClusterBytes;
};

} // namespace ybi
