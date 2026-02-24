#pragma once

#include "bvh/usd_subdiv_select.h"
#include "tesselation/subdivision_patch_types.h"

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

EdgeMapChecks RunEdgeMapChecks(const SelectedSubdivMesh &m,
                               const std::vector<SubdivisionPatch> &patches,
                               const SubdivisionEdgeMap &edgeMap);
