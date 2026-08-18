#pragma once

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
