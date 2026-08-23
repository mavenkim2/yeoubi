#pragma once

#include "core/vec3f.h"
#include "scene/material.h"
#include <cstdint>

#define OPENPBR_PARAMETERS(X) \
    X(base_weight, float) \
    X(base_color, Vec3f) \
    X(base_metalness, float) \
    X(base_diffuse_roughness, float) \
    X(specular_weight, float) \
    X(specular_color, Vec3f) \
    X(specular_roughness, float) \
    X(specular_roughness_anisotropy, float) \
    X(specular_ior, float) \
    X(transmission_weight, float) \
    X(transmission_color, Vec3f) \
    X(transmission_depth, float) \
    X(transmission_scatter, Vec3f) \
    X(transmission_scatter_anisotropy, float) \
    X(transmission_dispersion_scale, float) \
    X(transmission_dispersion_abbe_number, float) \
    X(subsurface_weight, float) \
    X(subsurface_color, Vec3f) \
    X(subsurface_radius, float) \
    X(subsurface_radius_scale, Vec3f) \
    X(subsurface_scatter_anisotropy, float) \
    X(coat_weight, float) \
    X(coat_color, Vec3f) \
    X(coat_roughness, float) \
    X(coat_roughness_anisotropy, float) \
    X(coat_ior, float) \
    X(coat_darkening, float) \
    X(fuzz_weight, float) \
    X(fuzz_color, Vec3f) \
    X(fuzz_roughness, float) \
    X(emission_luminance, float) \
    X(emission_color, Vec3f) \
    X(thin_film_weight, float) \
    X(thin_film_thickness, float) \
    X(thin_film_ior, float) \
    X(geometry_opacity, float) \
    X(geometry_thin_walled, bool) \
    X(geometry_normal, Vec3f) \
    X(geometry_tangent, Vec3f) \
    X(geometry_coat_normal, Vec3f) \
    X(geometry_coat_tangent, Vec3f)

namespace ybi
{

struct OpenPbrMaterial
{
    MaterialGraph graph;

#define DECLARE(name, type) MaterialNodeId name;
    OPENPBR_PARAMETERS(DECLARE)
#undef DECLARE
};

} // namespace ybi
