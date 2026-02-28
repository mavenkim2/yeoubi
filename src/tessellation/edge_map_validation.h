#pragma once

#include "scene/subdivision_mesh.h"
#include "tessellation/subdivision_patch_types.h"

#include <vector>

struct EdgeMapChecks
{
    int missingNgonMidpoints = 0;
    int duplicateMidpointVertices = 0;
    int missingGeneratedPatchEdges = 0;
    int badBoundaryFlags = 0;
    int badNonManifoldFlags = 0;
    bool ok = true;
};

EdgeMapChecks RunEdgeMapChecks(const ybi::SubdivisionMesh &m,
                               const std::vector<SubdivisionPatch> &patches,
                               const SubdivisionEdgeMap &edgeMap);
