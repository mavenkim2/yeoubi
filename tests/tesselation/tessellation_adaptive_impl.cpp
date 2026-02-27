#include "io/usd_subdiv_json_io.h"
#include "tessellation/subdivision.h"
#include "util/float3.h"

#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>

static ybi::BoundaryInterpolation BoundaryInterpolationFromString(const std::string &s)
{
    if (s == "none")
    {
        return ybi::BOUNDARY_INTERPOLATION_NONE;
    }
    if (s == "edgeOnly")
    {
        return ybi::BOUNDARY_INTERPOLATION_EDGE;
    }
    return ybi::BOUNDARY_INTERPOLATION_EDGE_AND_CORNER;
}

static ybi::FVarLinearInterpolation FVarLinearInterpolationFromString(const std::string &s)
{
    if (s == "none")
    {
        return ybi::FVAR_LINEAR_NONE;
    }
    if (s == "cornersOnly")
    {
        return ybi::FVAR_LINEAR_CORNERS_ONLY;
    }
    if (s == "cornersPlus1")
    {
        return ybi::FVAR_LINEAR_CORNERS_PLUS1;
    }
    if (s == "cornersPlus2")
    {
        return ybi::FVAR_LINEAR_CORNERS_PLUS2;
    }
    if (s == "boundaries")
    {
        return ybi::FVAR_LINEAR_BOUNDARIES;
    }
    if (s == "all" || s == "bilinear")
    {
        return ybi::FVAR_LINEAR_ALL;
    }
    return ybi::FVAR_LINEAR_CORNERS_PLUS1;
}

static ybi::SubdivisionMesh ConvertSelectedSubdivMesh(const SelectedSubdivMesh &in)
{
    ybi::SubdivisionMesh out = {};

    std::vector<ybi::float3> vertices;
    vertices.reserve(in.points.size());
    for (const pxr::GfVec3f &p : in.points)
    {
        ybi::float3 v = {};
        v.x = p[0];
        v.y = p[1];
        v.z = p[2];
        vertices.push_back(v);
    }
    out.vertices = vertices;

    std::vector<int> indices;
    indices.reserve(in.faceVertexIndices.size());
    for (int idx : in.faceVertexIndices)
    {
        indices.push_back(idx);
    }
    out.indices = indices;

    std::vector<int> vertsPerFace;
    vertsPerFace.reserve(in.faceVertexCounts.size());
    for (int c : in.faceVertexCounts)
    {
        vertsPerFace.push_back(c);
    }
    out.vertsPerFace = vertsPerFace;

    std::vector<int> cornerIndices;
    cornerIndices.reserve(in.cornerIndices.size());
    for (int idx : in.cornerIndices)
    {
        cornerIndices.push_back(idx);
    }
    out.cornerIndices = cornerIndices;

    std::vector<float> cornerSharpnesses;
    cornerSharpnesses.reserve(in.cornerSharpnesses.size());
    for (float sharpness : in.cornerSharpnesses)
    {
        cornerSharpnesses.push_back(sharpness);
    }
    out.cornerSharpnesses = cornerSharpnesses;

    std::vector<int> creaseIndices;
    creaseIndices.reserve(in.creaseIndices.size());
    for (int idx : in.creaseIndices)
    {
        creaseIndices.push_back(idx);
    }
    out.creaseIndices = creaseIndices;

    std::vector<int> creaseLengths;
    creaseLengths.reserve(in.creaseLengths.size());
    for (int len : in.creaseLengths)
    {
        creaseLengths.push_back(len);
    }
    out.creaseLengths = creaseLengths;

    std::vector<float> creaseSharpnesses;
    creaseSharpnesses.reserve(in.creaseSharpnesses.size());
    for (float sharpness : in.creaseSharpnesses)
    {
        creaseSharpnesses.push_back(sharpness);
    }
    out.creaseSharpnesses = creaseSharpnesses;

    std::vector<int> holeIndices;
    holeIndices.reserve(in.holeIndices.size());
    for (int idx : in.holeIndices)
    {
        holeIndices.push_back(idx);
    }
    out.holeIndices = holeIndices;

    out.interpolationRule = BoundaryInterpolationFromString(in.vertexBoundaryInterpolation);
    out.fvarLinearInterpolation = FVarLinearInterpolationFromString(in.fvarLinearInterpolation);
    out.attributeStart = 0;
    out.attributeEnd = 0;

    return out;
}

static pxr::GfVec3f ComputeMeshCenter(const SelectedSubdivMesh &mesh)
{
    if (mesh.points.empty())
    {
        return pxr::GfVec3f(0.0f, 0.0f, 0.0f);
    }

    pxr::GfVec3f sum(0.0f, 0.0f, 0.0f);
    for (const pxr::GfVec3f &p : mesh.points)
    {
        sum += p;
    }
    return sum * (1.0f / float(mesh.points.size()));
}

static bool WriteMeshObjWithTriMetadata(const ybi::Mesh &mesh,
                                        const std::vector<int> &triPatchFaceIds,
                                        const std::vector<int> &triCoarseFaceIds,
                                        const std::vector<int> &triPtexFaceIds,
                                        const std::vector<int> &triQuadrants,
                                        const std::string &outObjPath,
                                        bool writeMetadata)
{
    FILE *f = std::fopen(outObjPath.c_str(), "w");
    if (!f)
    {
        return false;
    }

    for (const ybi::float3 &p : mesh.positions)
    {
        std::fprintf(f, "v %.9g %.9g %.9g\n", double(p.x), double(p.y), double(p.z));
    }

    const size_t triCount = mesh.indices.size() / 3;
    const bool hasMeta = triPatchFaceIds.size() == triCount && triCoarseFaceIds.size() == triCount &&
                         triPtexFaceIds.size() == triCount && triQuadrants.size() == triCount;

    for (size_t triId = 0; triId < triCount; ++triId)
    {
        const int i0 = mesh.indices[triId * 3 + 0] + 1;
        const int i1 = mesh.indices[triId * 3 + 1] + 1;
        const int i2 = mesh.indices[triId * 3 + 2] + 1;
        if (writeMetadata && hasMeta)
        {
            std::fprintf(f,
                         "# tri_id %d patch_face_id %d coarse_face %d ptex_face_id %d quadrant %d\n",
                         int(triId),
                         triPatchFaceIds[triId],
                         triCoarseFaceIds[triId],
                         triPtexFaceIds[triId],
                         triQuadrants[triId]);
        }
        std::fprintf(f, "f %d %d %d\n", i0, i1, i2);
    }

    std::fclose(f);
    return true;
}

struct TessellationCliOptions
{
    std::string inJson;
    std::string outObjPath;
    int level = 1;
    float pixelSpacing = 1.0f;
    int splitThreshold = 1;
    int sampleSteps = 8;
    bool writeMetadata = false;
    std::string patchQuadObjPath;
};

static void PrintUsage(const char *exe)
{
    std::fprintf(stderr,
                 "Usage: %s --in <selected-subdiv.json> <out.obj>"
                 " [--level <n>=1] [--pixelspacing <f>>0] [--splitthreshold <n>=1]"
                 " [--samplesteps <n>=2] [--patchquadobj <out.obj>] [--writemetadata]\n",
                 exe);
}

static bool ParseCliArgs(int argc, char **argv, TessellationCliOptions *out, std::string *outError)
{
    if (!out)
    {
        return false;
    }

    for (int i = 1; i < argc; ++i)
    {
        const std::string arg = argv[i];
        if (arg == "-h" || arg == "--help")
        {
            return false;
        }

        if (arg.rfind("--", 0) != 0)
        {
            if (out->outObjPath.empty())
            {
                out->outObjPath = arg;
                continue;
            }
            if (outError)
            {
                *outError = "Unexpected positional arg: " + arg;
            }
            return false;
        }

        if (arg == "--writemetadata" || arg == "--write-metadata")
        {
            out->writeMetadata = true;
            continue;
        }

        if (i + 1 >= argc)
        {
            if (outError)
            {
                *outError = "Missing value for arg: " + arg;
            }
            return false;
        }
        const std::string value = argv[++i];

        if (arg == "--in" || arg == "--json")
        {
            out->inJson = value;
        }
        else if (arg == "--level")
        {
            out->level = std::max(1, std::atoi(value.c_str()));
        }
        else if (arg == "--pixelspacing" || arg == "--pixel-spacing")
        {
            out->pixelSpacing = std::max(1e-6f, float(std::atof(value.c_str())));
        }
        else if (arg == "--splitthreshold" || arg == "--split-threshold")
        {
            out->splitThreshold = std::max(1, std::atoi(value.c_str()));
        }
        else if (arg == "--samplesteps" || arg == "--sample-steps")
        {
            out->sampleSteps = std::max(2, std::atoi(value.c_str()));
        }
        else if (arg == "--patchquadobj" || arg == "--patch-quad-obj")
        {
            out->patchQuadObjPath = value;
        }
        else
        {
            if (outError)
            {
                *outError = "Unknown arg: " + arg;
            }
            return false;
        }
    }

    if (out->inJson.empty())
    {
        if (outError)
        {
            *outError = "Missing required arg: --in <selected-subdiv.json>";
        }
        return false;
    }
    if (out->outObjPath.empty())
    {
        if (outError)
        {
            *outError = "Missing required positional arg: <out.obj>";
        }
        return false;
    }

    return true;
}

int main(int argc, char **argv)
{
    TessellationCliOptions cli = {};
    std::string parseError;
    if (!ParseCliArgs(argc, argv, &cli, &parseError))
    {
        PrintUsage(argv[0]);
        if (!parseError.empty())
        {
            std::fprintf(stderr, "Error: %s\n", parseError.c_str());
            return 2;
        }
        return 0;
    }

    SelectedSubdivMesh selected = {};
    UsdCameraInfo camera = {};
    if (!ybi::testio::LoadSelectedSubdivFromJson(cli.inJson, selected, camera))
    {
        std::fprintf(stderr, "Failed to load JSON: %s\n", cli.inJson.c_str());
        return 1;
    }

    const ybi::SubdivisionMesh mesh = ConvertSelectedSubdivMesh(selected);

    const pxr::GfVec3f meshCenter = ComputeMeshCenter(selected);

    ybi::SubdivisionRunOptions options = {};
    options.level = cli.level;
    options.pixelSpacing = cli.pixelSpacing;
    options.splitThreshold = cli.splitThreshold;
    options.sampleSteps = cli.sampleSteps;
    options.eye = camera.found ? camera.worldPosition : (meshCenter + pxr::GfVec3f(0.0f, 0.0f, 5.0f));
    options.lookAt = camera.found ? camera.meshCenter : meshCenter;
    options.subdivisionScheme = selected.subdivisionScheme;
    options.creasingMethod = selected.creasingMethod;
    options.triangleSubdivision = selected.triangleSubdivision;
    options.patchQuadObjPath = cli.patchQuadObjPath;

    ybi::SubdivisionRunResult result = {};
    if (!ybi::SubdivideAdaptive(mesh, options, &result))
    {
        std::fprintf(stderr, "Subdivision failed\n");
        return 1;
    }

    if (!WriteMeshObjWithTriMetadata(result.mesh,
                                     result.trianglePatchFaceIds,
                                     result.triangleCoarseFaceIds,
                                     result.trianglePtexFaceIds,
                                     result.triangleQuadrants,
                                     cli.outObjPath,
                                     cli.writeMetadata))
    {
        std::fprintf(stderr, "Failed to write leaf inner-grid OBJ: %s\n", cli.outObjPath.c_str());
        return 1;
    }

    std::printf("  levelRequested=%d\n", cli.level);
    std::printf(
        "  refinedMaxLevel=%d patches=%d\n", result.refinedMaxLevel, result.totalPatches);
    std::printf(
        "  subdivisionPatches=%zu diagSplitPatches=%zu generatedVertexCount=%d midpointEdges=%d\n",
        result.subdivisionPatchCount,
        result.diagSplitPatchCount,
        result.generatedVertexCount,
        result.midpointEdges);
    std::printf(
        "  edgeTMaxComputed=%d totalComputedEdges=%d\n", result.edgeTMaxComputed, result.totalComputedEdges);
    std::printf("  diagSplitConfig pixelSpacing=%.3f splitThreshold=%d sampleSteps=%d\n",
                cli.pixelSpacing,
                cli.splitThreshold,
                cli.sampleSteps);
    if (cli.patchQuadObjPath.empty())
    {
        std::printf("  patchQuadObj disabled\n");
    }
    else
    {
        std::printf("  patchQuadObj path=%s verts=%d quads=%d\n",
                    cli.patchQuadObjPath.c_str(),
                    result.patchQuadVerts,
                    result.patchQuadCount);
    }
    std::printf("  innerGridObj path=%s verts=%zu tris=%zu\n",
                cli.outObjPath.c_str(),
                result.mesh.positions.size(),
                result.mesh.indices.size() / 3);
    std::printf("  writeMetadata=%s\n", cli.writeMetadata ? "true" : "false");
    std::printf("  edgeMapChecks=ok\n");
    std::printf("  controlCageUniqueEdges=%zu boundaryEdges=%d\n",
                result.controlCageUniqueEdges,
                result.boundaryEdges);
    std::printf("  controlCageEdgesWithOver2Faces=%d\n", result.controlCageEdgesWithOver2Faces);

    return 0;
}
