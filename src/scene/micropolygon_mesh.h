#pragma once

#include "scene/subdivision.h"
#include "util/vec3.h"

namespace ybi
{

struct Grid
{
    int gridIndexStart;
    int width;
    int height;
};

struct MicropolygonMesh
{
    Array<Vec3> positions;
    Array<Grid> grids;
};

} // namespace ybi
