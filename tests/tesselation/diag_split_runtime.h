#pragma once

#include "io/usd_subdiv_json_io.h"
#include "tesselation/edge_rate_obj_io.h"

#include <opensubdiv/far/patchMap.h>
#include <opensubdiv/far/patchTable.h>
#include <opensubdiv/far/ptexIndices.h>
#include <opensubdiv/far/topologyRefiner.h>

#include <array>
#include <vector>

namespace ybi::tesselation
{

struct Primvar3
{
    float x = 0.0f;
    float y = 0.0f;
    float z = 0.0f;
    Primvar3() = default;
    explicit Primvar3(const pxr::GfVec3f &p);
    void Clear();
    void AddWithWeight(const Primvar3 &p, float w);
    pxr::GfVec3f ToGf() const;
};

struct PtexFaceAdj
{
    std::array<int, 4> adjFace = {-1, -1, -1, -1};
    std::array<int, 4> adjEdge = {-1, -1, -1, -1};
    std::array<int, 4> edgeIndex = {-1, -1, -1, -1};
    bool fromNgon = false;
};

struct EdgeFactorCamera
{
    pxr::GfVec3f origin = pxr::GfVec3f(0.0f);
    pxr::GfVec3f target = pxr::GfVec3f(0.0f);
    pxr::GfVec3f up = pxr::GfVec3f(0.0f, 1.0f, 0.0f);
    float fovYDegrees = 50.0f;
    int width = 1920;
    int height = 1080;
};

struct EdgeFactorSettings
{
    int numSamples = 4;
    float pixelSpacing = 1.0f;
    int minRate = 1;
    int maxRate = 32;
};

struct RefinedPositions
{
    std::vector<Primvar3> values;
    std::vector<int> levelStarts;
    int numRefinerVerts = 0;
};

struct TriMesh
{
    std::vector<pxr::GfVec3f> positions;
    std::vector<int> indices;
};

constexpr int kDiagSplitNonUniform = -1;

struct DiagEdgeState
{
    int child0 = -1;
    int child1 = -1;
    int rate = kDiagSplitNonUniform;
    float firstT0 = 0.0f;
    float firstT1 = 0.0f;
    int firstFace = -1;
    bool firstSet = false;
};

struct DiagSubPatch
{
    int ptexFace = -1;
    std::array<int, 4> edgeIndex = {-1, -1, -1, -1};
    float u0 = 0.0f;
    float v0 = 0.0f;
    float u1 = 1.0f;
    float v1 = 1.0f;
    int depth = 0;
};

struct DiagSplitBuildResult
{
    std::vector<DiagSubPatch> patches;
    int splitCount = 0;
    int maxDepthReached = 0;
    int skippedPtexFaces = 0;
};

EdgeFactorCamera BuildEdgeFactorCamera(const pxr::VtVec3fArray &points,
                                       const UsdCameraInfo &cameraInfo,
                                       float cameraDistanceScale);
RefinedPositions BuildRefinedPositions(const OpenSubdiv::Far::TopologyRefiner &refiner,
                                       const OpenSubdiv::Far::PatchTable &patchTable,
                                       const pxr::VtVec3fArray &points);
int BuildUniquePtexEdgeIds(const OpenSubdiv::Far::TopologyRefiner &refiner,
                           const OpenSubdiv::Far::TopologyLevel &level0,
                           const OpenSubdiv::Far::PtexIndices &ptex,
                           std::vector<PtexFaceAdj> &facesOut,
                           int &numPtexFacesOut);
std::vector<DiagEdgeState> ComputeEdgeFactors(const OpenSubdiv::Far::PatchMap &patchMap,
                                              const OpenSubdiv::Far::PatchTable &patchTable,
                                              const std::vector<Primvar3> &positions,
                                              const std::vector<PtexFaceAdj> &faces,
                                              int numPtexFaces,
                                              int uniqueEdgeCount,
                                              const EdgeFactorCamera &camera,
                                              const EdgeFactorSettings &settings);
void ApplyNgonNonUniformConstraint(const std::vector<PtexFaceAdj> &faces,
                                   int numPtexFaces,
                                   std::vector<DiagEdgeState> &edgeRates);
std::vector<ybi::testio::EdgeRateDebugLine> BuildEdgeRateDebugLines(
    const OpenSubdiv::Far::PatchMap &patchMap,
    const OpenSubdiv::Far::PatchTable &patchTable,
    const std::vector<Primvar3> &positions,
    const std::vector<PtexFaceAdj> &faces,
    int numPtexFaces,
    const std::vector<DiagEdgeState> &edgeFactors);
DiagSplitBuildResult BuildDiagSplitSubPatches(const OpenSubdiv::Far::PatchMap &patchMap,
                                              const OpenSubdiv::Far::PatchTable &patchTable,
                                              const std::vector<Primvar3> &positions,
                                              const std::vector<PtexFaceAdj> &faces,
                                              int numPtexFaces,
                                              const EdgeFactorCamera &camera,
                                              const EdgeFactorSettings &settings,
                                              std::vector<DiagEdgeState> &edges);
TriMesh TessellateDiagSplitSubPatches(const OpenSubdiv::Far::PatchMap &patchMap,
                                      const OpenSubdiv::Far::PatchTable &patchTable,
                                      const std::vector<Primvar3> &positions,
                                      const std::vector<DiagSubPatch> &patches,
                                      const std::vector<DiagEdgeState> &edges,
                                      int fallbackEdgeRate,
                                      int &skippedSubPatchesOut);

} // namespace ybi::tesselation
