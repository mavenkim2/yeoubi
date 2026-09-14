#pragma once

#include <cstdint>

namespace ybi
{

#define RENDERER_NODE_DEFINITIONS(X) \
    X(ConstantFloat, Float) \
    X(ConstantBool, Bool) \
    X(ConstantColor, Vec3) \
    X(ImageVector3, Vec3) \
    X(ImageFloat, Float) \
    X(AttributeValueVector2, Vec2) \
    X(HeightToNormalVector3, Vec3) \
    X(InvertFloat, Float) \
    X(ExtractVector3, Float) \
    X(Combine3Vector3, Vec3) \
    X(NormalMapVector2, Vec3) \
    X(ColorCorrectColor3, Vec3) \
    X(RemapColor3, Vec3) \
    X(Ramp4Float, Float) \
    X(Separate3Color3, Float)

#define MATERIALX_DEFINITIONS(X) \
    X(constant_boolean, ConstantBool) \
    X(image_color3, ImageVector3) \
    X(image_vector3, ImageVector3) \
    X(image_float, ImageFloat) \
    X(geompropvalue_vector2, AttributeValueVector2) \
    X(heighttonormal_vector3, HeightToNormalVector3) \
    X(invert_float, InvertFloat) \
    X(extract_vector3, ExtractVector3) \
    X(extract_color3, ExtractVector3) \
    X(combine3_vector3, Combine3Vector3) \
    X(normalmap_vector2, NormalMapVector2) \
    X(colorcorrect_color3, ColorCorrectColor3) \
    X(remap_color3, RemapColor3) \
    X(ramp4_float, Ramp4Float) \
    X(separate3_color3, Separate3Color3)

enum class MaterialValueType : uint8_t
{
    None,
    Bool,
    Int,
    Float,
    Vec2,
    Vec3,
    Vec4,
};

constexpr uint8_t GetMaterialValueRegisterCount(MaterialValueType type)
{
    switch (type)
    {
        case MaterialValueType::Bool:
        case MaterialValueType::Int:
        case MaterialValueType::Float:
            return 1;
        case MaterialValueType::Vec2:
            return 2;
        case MaterialValueType::Vec3:
            return 3;
        case MaterialValueType::Vec4:
            return 4;
        case MaterialValueType::None:
            return 0;
    }

    return 0;
}

enum class MaterialOpcode : uint8_t
{
    Invalid = 0,

#define ENUM(name, outputType) name,
    RENDERER_NODE_DEFINITIONS(ENUM)
#undef ENUM

        Count,
};

constexpr const char *ConvertMaterialOpcodeToString(MaterialOpcode opcode)
{
    switch (opcode)
    {
#define CASE(name, outputType) \
    case MaterialOpcode::name: \
        return #name;
        RENDERER_NODE_DEFINITIONS(CASE)
#undef CASE

        case MaterialOpcode::Invalid:
            return "Invalid";
        case MaterialOpcode::Count:
            return "Count";
    }

    return "Invalid";
}

constexpr MaterialValueType GetMaterialOpcodeOutputType(MaterialOpcode opcode)
{
    switch (opcode)
    {
#define CASE(name, outputType) \
    case MaterialOpcode::name: \
        return MaterialValueType::outputType;
        RENDERER_NODE_DEFINITIONS(CASE)
#undef CASE

        case MaterialOpcode::Invalid:
        case MaterialOpcode::Count:
            return MaterialValueType::None;
    }

    return MaterialValueType::None;
}

constexpr uint8_t GetMaterialOpcodeOutputRegisterCount(MaterialOpcode opcode)
{
    return GetMaterialValueRegisterCount(GetMaterialOpcodeOutputType(opcode));
}

constexpr bool IsTextureOpcode(MaterialOpcode opcode)
{
    return opcode == MaterialOpcode::ImageVector3 || opcode == MaterialOpcode::ImageFloat;
}

constexpr bool IsAttributeOpcode(MaterialOpcode opcode)
{
    return opcode == MaterialOpcode::AttributeValueVector2;
}

} // namespace ybi
