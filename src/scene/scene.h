#pragma once

#include <vector>

YBI_NAMESPACE_BEGIN

struct Mesh
{
    float3 *positions;
    int *indices;
};

struct Scene
{
    std::vector<Mesh> meshes;
};

YBI_NAMESPACE_END
