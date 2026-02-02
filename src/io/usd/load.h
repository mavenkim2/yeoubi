#pragma once

#include <pxr/usd/usd/stage.h>
#include <pxr/usd/usdGeom/mesh.h>

YBI_NAMESPACE_BEGIN

static void Test() { pxr::UsdStageRefPtr stage = pxr::UsdStage::Open(filePath); }

YBI_NAMESPACE_END
