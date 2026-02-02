#pragma once

#include <pxr/usd/usd/stage.h>
#include <pxr/usd/usdGeom/mesh.h>
#include <vector>

YBI_NAMESPACE_BEGIN

template <typename T, typename V>
using vector = std::vector<T, V>;

static void Test(std::string filePath)
{
    pxr::UsdStageRefPtr stage = pxr::UsdStage::Open(filePath);
}

YBI_NAMESPACE_END
