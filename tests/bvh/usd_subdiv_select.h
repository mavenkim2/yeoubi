#pragma once

#include <pxr/usd/usd/stage.h>
#include <pxr/base/vt/array.h>
#include <pxr/usd/sdf/path.h>
#include <string>

struct SelectedSubdivMesh
{
    pxr::SdfPath path;
    pxr::VtVec3fArray points;
    pxr::VtIntArray faceVertexCounts;
    pxr::VtIntArray faceVertexIndices;
    std::string subdivisionScheme;
};

bool SelectLargestCatmullClarkMesh(const pxr::UsdStageRefPtr &stage, SelectedSubdivMesh &meshOut);
