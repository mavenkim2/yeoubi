#pragma once

#include <cstdint>

namespace ybi
{

#define RENDERER_NODE_DEFINITIONS(X) \
    X(ConstantFloat, Float) \
    X(ConstantBool, Bool) \
    X(ConstantColor, Vec3) \
    X(ConstantVector2, Vec2) \
    X(ImageVector2, Vec2) \
    X(ImageVector3, Vec3) \
    X(ImageVector4, Vec4) \
    X(ImageFloat, Float) \
    X(AttributeValueVector2, Vec2) \
    X(HeightToNormalVector3, Vec3) \
    X(InvertFloat, Float) \
    X(ExtractVector3, Float) \
    X(Combine3Vector3, Vec3) \
    X(NormalMapVector2, Vec3) \
    X(ColorCorrectColor3, Vec3) \
    X(RemapVector3, Vec3) \
    X(Ramp4Float, Float) \
    X(Separate3Color3, Float) \
    X(Mix, Float)

// Value nodes compiled into MaterialX callables rather than the legacy shader bytecode.
#define MATERIALX_VALUE_NODE_DEFINITIONS(X) \
    X(TiledImageVector2, Vec2) \
    X(TiledImageVector3, Vec3) \
    X(TiledImageVector4, Vec4) \
    X(TiledImageFloat, Float) \
    X(ViewDirectionVector3, Vec3) \
    X(NormalVector3, Vec3) \
    X(DotProductVector3, Float) \
    X(AbsValFloat, Float) \
    X(AsinFloat, Float) \
    X(DivideFloat, Float) \
    X(DivideVector3, Vec3) \
    X(MultiplyVector3, Vec3) \
    X(ExpVector3, Vec3) \
    X(LnVector3, Vec3) \
    X(ConvertFloatVector2, Vec2) \
    X(ConvertFloatVector3, Vec3) \
    X(ConvertVector3Color3, Vec3) \
    X(ConvertColor3Vector3, Vec3) \
    X(RemapFloat, Float) \
    X(AddFloat, Float) \
    X(AddVector3, Vec3) \
    X(NormalizeVector3, Vec3) \
    X(MixFloat, Float) \
    X(ContrastFloat, Float) \
    X(ClampFloat, Float) \
    X(NormalMapFloat, Vec3) \
    X(ArtisticIorIor, Vec3) \
    X(ArtisticIorExtinction, Vec3) \
    X(PowerFloat, Float) \
    X(MixVector3, Vec3) \
    X(MixVector2, Vec2) \
    X(MultiplyFloat, Float)

#define MATERIALX_DEFINITIONS(X) \
    X(constant_boolean, ConstantBool) \
    X(constant_float, ConstantFloat) \
    X(constant_vector2, ConstantVector2) \
    X(image_vector2, ImageVector2) \
    X(image_color3, ImageVector3) \
    X(image_vector3, ImageVector3) \
    X(image_color4, ImageVector4) \
    X(image_vector4, ImageVector4) \
    X(image_float, ImageFloat) \
    X(tiledimage_vector2, TiledImageVector2) \
    X(tiledimage_vector3, TiledImageVector3) \
    X(tiledimage_color3, TiledImageVector3) \
    X(tiledimage_vector4, TiledImageVector4) \
    X(tiledimage_color4, TiledImageVector4) \
    X(tiledimage_float, TiledImageFloat) \
    X(geompropvalue_vector2, AttributeValueVector2) \
    X(heighttonormal_vector3, HeightToNormalVector3) \
    X(invert_float, InvertFloat) \
    X(extract_vector3, ExtractVector3) \
    X(extract_color3, ExtractVector3) \
    X(combine3_vector3, Combine3Vector3) \
    X(normalmap_vector2, NormalMapVector2) \
    X(normalmap, NormalMapFloat) \
    X(normalmap_float, NormalMapFloat) \
    X(colorcorrect_color3, ColorCorrectColor3) \
    X(remap_color3, RemapVector3) \
    X(remap_vector3, RemapVector3) \
    X(ramp4_float, Ramp4Float) \
    X(separate3_color3, Separate3Color3) \
    X(mix_float, MixFloat) \
    X(mix_vector3, MixVector3) \
    X(mix_vector2, MixVector2) \
    X(constant_color3, ConstantColor) \
    X(constant_vector3, ConstantColor) \
    X(viewdirection_vector3, ViewDirectionVector3) \
    X(normal_vector3, NormalVector3) \
    X(dotproduct_vector3, DotProductVector3) \
    X(absval_float, AbsValFloat) \
    X(asin_float, AsinFloat) \
    X(power_float, PowerFloat) \
    X(divide_float, DivideFloat) \
    X(divide_vector3, DivideVector3) \
    X(multiply_vector3, MultiplyVector3) \
    X(multiply_color3, MultiplyVector3) \
    X(multiply_float, MultiplyFloat) \
    X(exp_vector3, ExpVector3) \
    X(ln_vector3, LnVector3) \
    X(convert_float_vector2, ConvertFloatVector2) \
    X(convert_float_vector3, ConvertFloatVector3) \
    X(convert_vector3_color3, ConvertVector3Color3) \
    X(convert_color3_vector3, ConvertColor3Vector3) \
    X(remap_float, RemapFloat) \
    X(add_float, AddFloat) \
    X(add_vector3, AddVector3) \
    X(normalize_vector3, NormalizeVector3) \
    X(contrast_float, ContrastFloat) \
    X(clamp_float, ClampFloat)

#define MATERIALX_BSDFS(X) \
    X(layer) \
    X(mix) \
    X(sheen) \
    X(oren_nayar_diffuse) \
    X(dielectric) \
    X(conductor)

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

constexpr const char *GetMaterialValueTypeString(MaterialValueType type)
{
    switch (type)
    {
        case MaterialValueType::Bool:
        {
            return "bool";
        }
        case MaterialValueType::Int:
        {
            return "int";
        }
        case MaterialValueType::Float:
        {
            return "float";
        }
        case MaterialValueType::Vec2:
        {
            return "Vec2f";
        }
        case MaterialValueType::Vec3:
        {
            return "Vec3f";
        }
        case MaterialValueType::Vec4:
        {
            return "Vec4f";
        }
        case MaterialValueType::None:
        {
            return "";
        }
    }
    return "";
}

enum class MaterialOpcode : uint8_t
{
    Invalid = 0,

#define ENUM(name, outputType) name,
    RENDERER_NODE_DEFINITIONS(ENUM) MATERIALX_VALUE_NODE_DEFINITIONS(ENUM)
#undef ENUM
#define ENUM(name) name##_bsdf,
        MATERIALX_BSDFS(ENUM)
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
        MATERIALX_VALUE_NODE_DEFINITIONS(CASE)
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
        MATERIALX_VALUE_NODE_DEFINITIONS(CASE)
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

constexpr const char *GetMaterialOpcodeOutputString(MaterialOpcode opcode)
{
    return GetMaterialValueTypeString(GetMaterialOpcodeOutputType(opcode));
}

constexpr bool IsTextureOpcode(MaterialOpcode opcode)
{
    return opcode == MaterialOpcode::ImageFloat || opcode == MaterialOpcode::ImageVector2 ||
           opcode == MaterialOpcode::ImageVector3 || opcode == MaterialOpcode::ImageVector4 ||
           opcode == MaterialOpcode::TiledImageFloat ||
           opcode == MaterialOpcode::TiledImageVector2 ||
           opcode == MaterialOpcode::TiledImageVector3 ||
           opcode == MaterialOpcode::TiledImageVector4;
}

constexpr bool IsAttributeOpcode(MaterialOpcode opcode)
{
    return opcode == MaterialOpcode::AttributeValueVector2;
}

constexpr bool IsConstantOpcode(MaterialOpcode opcode)
{
    return opcode == MaterialOpcode::ConstantFloat || opcode == MaterialOpcode::ConstantBool ||
           opcode == MaterialOpcode::ConstantColor || opcode == MaterialOpcode::ConstantVector2;
}

} // namespace ybi
