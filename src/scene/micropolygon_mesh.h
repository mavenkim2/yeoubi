#pragma once

#include "scene/subdivision.h"
#include "util/float3.h"

YBI_NAMESPACE_BEGIN

struct Grid
{
    int gridIndexStart;
    int width;
    int height;
};

struct MicropolygonMesh
{
    Array<float3> positions;
    Array<Grid> grids;
};

YBI_NAMESPACE_END
