#pragma once

#include <pxr/pxr.h>
#include <pxr/usd/usdShade/shader.h>

PXR_NAMESPACE_OPEN_SCOPE
class UsdShadeMaterial;
PXR_NAMESPACE_CLOSE_SCOPE

namespace ybi
{

pxr::UsdShadeShader GetMaterialSurfaceShader(const pxr::UsdShadeMaterial &material);
void ProcessUsdShader(const pxr::UsdShadeShader &shader, pxr::UsdTimeCode timeCode);

} // namespace ybi
