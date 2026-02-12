#pragma once

#include "scene/subdivision_mesh.h"
#include <opensubdiv/far/patchTable.h>
#include <opensubdiv/far/topologyRefinerFactory.h>

YBI_NAMESPACE_BEGIN

struct Scene;

struct SubdivisionResult
{
    std::shared_ptr<OpenSubdiv::Far::TopologyRefiner> refiner;
    std::unique_ptr<const OpenSubdiv::Far::PatchTable> patchTable;

    Array<float3> limitPositions;
};

struct MicropolygonMesh
{
    Array<float3> positions;
    Array<int> indices;
};

void Subdivision(Scene *scene, const SubdivisionMesh &mesh, int refineLevel = 1);
MicropolygonMesh CreateMicropolygons(const SubdivisionResult &subdivisionResult,
                                     const SubdivisionMesh &controlMesh,
                                     int edgeRateU,
                                     int edgeRateV);

YBI_NAMESPACE_END
