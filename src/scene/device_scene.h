#pragma once

#include "core/affine_transform.h"
#include "core/vec2f.h"
#include "core/vec3f.h"
#include "render/material_opcodes.h"
#include "scene/device_cluster.h"
#include "scene/openpbr_parameter_list.h"
#include <cstdint>

namespace ybi
{

using AttributeId = uint16_t;

enum class AttributeType
{
    None,
    Vertex,
    FaceVarying,
};

enum class MaterialHandleType : uint32_t
{
    Preview,
    OpenPbr,
    MaterialX,
};

struct MaterialHandle
{
    MaterialHandleType type = MaterialHandleType::Preview;
    uint32_t index = 0;
};

struct DeviceMaterial
{
    // NOTE: = ~0u if no udims
    uint32_t udimBase;
    uint32_t udimTableOffset;

    bool valid;
};

struct DeviceMaterialExpression
{
    uint32_t bytecodeOffset;
    uint32_t outputRegister;
};

struct DeviceMaterial2
{
#define DECLARE(name, type) DeviceMaterialExpression name;
    OPENPBR_PARAMETERS(DECLARE)
#undef DECLARE
};

struct DeviceMaterial3
{
    uint32_t sampleCallableIndex;
    uint32_t evaluateCallableIndex;
    uint32_t pdfCallableIndex;
};

struct DeviceAttribute
{
    AttributeId id;
    AttributeType attributeType;

    // NOTE: these are unused when YBI_USE_CLUSTERS is enabled
    int index;
    int indicesIndex;
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
    uint32_t primvarOffset;
    uint32_t primvarCount;

    int materialIndex = -1;
    int lightIndex = -1;
};

enum class LightType
{
    Rect,
};

struct DeviceLight
{
    // Rect
    AffineTransform worldFromLocal;
    AffineTransform localFromWorld;
    float width;
    float height;
    Vec3f radiance;
};

struct MaterialDescriptor
{
    // Indices are local to a child-before-parent descriptor array. For layer,
    // child0 is top and child1 is base; for mix, they are bg and fg.
    MaterialOpcode materialType = MaterialOpcode::Invalid;
    uint32_t child0 = 0xffffffffu;
    uint32_t child1 = 0xffffffffu;
    uint32_t parameterOffset = 0xffffffffu;
    // Only the lion's default BSDF modes are currently supported.
    uint32_t flags = 0;
};

struct DeviceScene
{
    Object *objects;
    DeviceMaterial *materials;
    DeviceMaterial2 *deviceMaterial2s;
    DeviceMaterial3 *deviceMaterial3s;
    MaterialHandle *materialHandles;

    Vec3f *positions;
    uint32_t *normals;
    Vec2f *uvs;
    int *uvIndices;
    int *triangleIndices;
    int *triangleMaterialIndices;

    DeviceCluster *compressedClusters;
    uint8_t *compressedClusterBytes;

    DeviceAttribute *deviceAttributes;
    uint32_t *shaderNodes;

    MaterialDescriptor *materialDescriptors;

    DeviceLight *deviceLights;
    uint32_t numLights;

    // TODO: populate these and use them when clusterizer isn't used
    int *attributeIndices;
    Vec2f *attributesVec2f;
};

} // namespace ybi
