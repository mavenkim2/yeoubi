#pragma once

#include "scene/openpbr_parameters.h"
#include <pxr/base/tf/token.h>
#include <pxr/pxr.h>
#include <pxr/usd/usdShade/shader.h>

PXR_NAMESPACE_OPEN_SCOPE
class UsdShadeMaterial;
PXR_NAMESPACE_CLOSE_SCOPE

namespace ybi
{

pxr::UsdShadeShader GetMaterialSurfaceShader(const pxr::UsdShadeMaterial &material);
void LoadOpenPbrMaterial(const pxr::UsdShadeShader &shader,
                         pxr::UsdTimeCode timeCode,
                         std::vector<std::string> &geomProps,
                         OpenPbrMaterial *material);

} // namespace ybi
