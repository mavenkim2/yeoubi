#pragma once

#include "core/matrix4x4.h"
#include "core/vec2f.h"
#include "core/vec3f.h"
#include "core/vec4f.h"
#include <string>
#include <variant>
#include <vector>

// TODO: instead of hardcoding,
#define MATERIALX_DEFINITIONS(X) \
    X(image_color3, Vec3f) \
    X(image_vector3, Vec3f) \
    X(image_float, float) \
    X(geompropvalue_vector2, Vec2f) \
    X(heighttonormal_vector3, Vec3f) \
    X(invert_float, float) \
    X(extract_vector3, float) \
    X(extract_color3, float) \
    X(combine3_vector3, Vec3f) \
    X(normalmap_vector2, Vec3f) \
    X(colorcorrect_color3, Vec3f) \
    X(remap_color3, Vec3f) \
    X(ramp4_float, float) \
    X(separate3_color3, float)

namespace ybi
{
enum class MtlxNodeType
{
    None,
#define ENUM(name, type) name,
    MATERIALX_DEFINITIONS(ENUM)
#undef ENUM
};

struct MaterialNode;

using MaterialInputValue = std::variant<std::shared_ptr<MaterialNode>,
                                        bool,
                                        int,
                                        float,
                                        Vec2f,
                                        Vec3f,
                                        Vec4f,
                                        Matrix4x4,
                                        std::string>;

struct MaterialInput
{
    std::string name;
    MaterialInputValue value;
};

struct MaterialNode
{
    MtlxNodeType type;
    std::vector<MaterialInput> inputs;
};

inline const char *ConvertMtlxNodeTypeToString(MtlxNodeType nodeType)
{
#define SWITCH(name, type) \
    if (nodeType == MtlxNodeType::name) \
        return #name;

    MATERIALX_DEFINITIONS(SWITCH)
    return "None";
#undef SWITCH
}

inline size_t GetMtlxNodeOutputSize(MtlxNodeType nodeType)
{
#define SWITCH(name, type) \
    if (nodeType == MtlxNodeType::name) \
        return sizeof(type);
    MATERIALX_DEFINITIONS(SWITCH)
    return 0;
#undef SWITCH
}

} // namespace ybi
