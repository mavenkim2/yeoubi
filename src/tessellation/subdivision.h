#pragma once

#include "scene/scene.h"
#include "scene/subdivision_mesh.h"
#include "util/vec3.h"

#include <string>
#include <vector>

namespace ybi
{

struct SubdivisionRunOptions
{
    int level = 3;
    int maxDiagSplitDepth = 12;
    float pixelSpacing = 1.0f;
    int splitThreshold = 1;
    int sampleSteps = 8;

    int viewportWidth = 1920;
    int viewportHeight = 1080;
    Float4x4 cameraFromWorld = Float4x4::Identity();
    Float4x4 clipFromCamera = Float4x4::Identity();

    std::string subdivisionScheme = "catmullClark";
    SubdivisionCreasingMethod creasingMethod = SUBDIVISION_CREASING_UNIFORM;
    bool generateTriangleMetadata = true;

    std::string patchQuadObjPath;
};

struct SubdivisionRunResult
{
    Mesh mesh;

    std::vector<int> trianglePatchFaceIds;
    std::vector<int> triangleCoarseFaceIds;
    std::vector<int> trianglePtexFaceIds;
    std::vector<int> triangleQuadrants;

    int refinedMaxLevel = 0;
    int totalPatches = 0;

    size_t subdivisionPatchCount = 0;
    size_t diagSplitPatchCount = 0;
    int generatedVertexCount = 0;
    int midpointEdges = 0;

    int edgeTMaxComputed = 0;
    int totalComputedEdges = 0;

    int patchQuadVerts = 0;
    int patchQuadCount = 0;
    int innerGridTriangleCount = 0;
    int stitchingTriangleCount = 0;

    size_t controlCageUniqueEdges = 0;
    int boundaryEdges = 0;
    int controlCageEdgesWithOver2Faces = 0;
};

bool SubdivideAdaptive(const SubdivisionMesh &mesh,
                      const SubdivisionRunOptions &options,
                      SubdivisionRunResult *outResult);

} // namespace ybi
