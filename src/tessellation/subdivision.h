#pragma once

#include "scene/scene.h"
#include "scene/subdivision_mesh.h"
#include "util/float3.h"

#include <string>
#include <vector>

YBI_NAMESPACE_BEGIN

struct SubdivisionRunOptions
{
    int level = 3;
    int maxDiagSplitDepth = 12;
    float pixelSpacing = 1.0f;
    int splitThreshold = 1;
    int sampleSteps = 8;

    float3 eye = make_float3(0.0f, 0.0f, 5.0f);
    float3 lookAt = make_float3(0.0f, 0.0f, 0.0f);
    int viewportWidth = 1920;
    int viewportHeight = 1080;
    float verticalFovDegrees = 45.0f;

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

YBI_NAMESPACE_END
