#include "scene.h"
#include "../util/float3.h"

YBI_NAMESPACE_BEGIN

Mesh::Mesh(float3 *positions, int *indices) : positions(positions), indices(indices) {}

YBI_NAMESPACE_END
