#pragma once

#include "device/kernel_params.h"
#include "render/shaders/material_shaders.h"
#include "render/shaders/shading.h"
#include "scene/device_scene.h"

namespace ybi
{

struct OpenPbrBSDF
{
    float ior;
    // float
};

YBI_CPU_GPU void InitializeOpenPbrBSDF(const KernelParams &params,
                                       const DeviceMaterial2 &material,
                                       const ShadingInfo &info,
                                       OpenPbrBSDF &bsdf)
{
    float specularIor = EvaluateShader<float>(params, info, material.specular_ior);
    float specularWeight = EvaluateShader<float>(params, info, material.specular_weight);
    Vec3f specularColor = EvaluateShader<Vec3f>(params, info, material.specular_color);
}

} // namespace ybi
