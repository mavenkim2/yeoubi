#pragma once

#include "scene/subdivision_mesh.h"
#include <opensubdiv/far/patchTable.h>
#include <opensubdiv/far/topologyRefinerFactory.h>

YBI_NAMESPACE_BEGIN

struct Scene;

void Subdivision(Scene *scene, const SubdivisionMesh &mesh, int refineLevel = 1);

YBI_NAMESPACE_END
