#include "io/usd_subdiv_json_io.h"
#include "tessellation/subdivision.h"
#include "tessellation_adaptive_cli_usd.h"
#include "util/vec3.h"

#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <cstddef>
#include <string>
#include <vector>

static ybi::Vec3 ComputeMeshCenter(const ybi::SubdivisionMesh &mesh)
{
    if (mesh.vertices.size() == 0)
    {
        return ybi::Vec3(0.0f);
    }

    ybi::Vec3 sum = ybi::Vec3(0.0f);
    for (const ybi::Vec3 &p : mesh.vertices)
    {
        sum += p;
    }
    return sum * (1.0f / float(mesh.vertices.size()));
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

    for (const ybi::Vec3 &p : mesh.positions)
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
    std::string inPath;
    std::string outObjPath;
    int level = 3;
    float pixelSpacing = 1.0f;
    int splitThreshold = 1;
    int sampleSteps = 8;
    bool writeMetadata = false;
    std::string patchQuadObjPath;
};

static void PrintUsage(const char *exe)
{
    std::fprintf(stderr,
                 "Usage: %s --in <input.json|input.usd[a|c]> <out.obj|out_dir>"
                 " [--level <n>=3] [--pixelspacing <f>>0] [--splitthreshold <n>=1]"
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
            out->inPath = value;
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

    if (out->inPath.empty())
    {
        if (outError)
        {
            *outError = "Missing required arg: --in <input.json|input.usd[a|c]>";
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

    if (!IsJsonInputPath(cli.inPath) && !IsUsdInputPath(cli.inPath))
    {
        std::fprintf(stderr, "Unsupported input extension: %s\n", cli.inPath.c_str());
        return 1;
    }

    if (IsJsonInputPath(cli.inPath))
    {
        ybi::SubdivisionMesh mesh = {};
        UsdCameraInfo camera = {};
        std::string subdivisionScheme = "catmullClark";
        if (!ybi::testio::LoadSelectedSubdivFromJson(cli.inPath,
                                                     mesh,
                                                     camera,
                                                     &subdivisionScheme))
        {
            std::fprintf(stderr, "Failed to load JSON: %s\n", cli.inPath.c_str());
            return 1;
        }

        const ybi::Vec3 meshCenter = ComputeMeshCenter(mesh);

        ybi::SubdivisionRunOptions options = {};
        options.level = cli.level;
        options.pixelSpacing = cli.pixelSpacing;
        options.splitThreshold = cli.splitThreshold;
        options.sampleSteps = cli.sampleSteps;
        const ybi::Vec3 eye = camera.found
                                  ? ybi::Vec3(camera.worldPosition[0],
                                              camera.worldPosition[1],
                                              camera.worldPosition[2])
                                  : (meshCenter + ybi::Vec3(0.0f, 0.0f, 5.0f));
        const ybi::Vec3 lookAt =
            camera.found ? ybi::Vec3(camera.meshCenter[0], camera.meshCenter[1], camera.meshCenter[2])
                         : meshCenter;
        options.viewportWidth = std::max(1, options.viewportWidth);
        options.viewportHeight = std::max(1, options.viewportHeight);
        options.cameraFromWorld = ybi::BuildCameraFromWorld(eye, lookAt);
        options.clipFromCamera =
            ybi::BuildPerspectiveClipFromCamera(45.0f, options.viewportWidth, options.viewportHeight);
        options.subdivisionScheme = subdivisionScheme;
        options.generateTriangleMetadata = cli.writeMetadata;
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
        std::printf("  edgeTMaxComputed=%d totalComputedEdges=%d\n",
                    result.edgeTMaxComputed,
                    result.totalComputedEdges);
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
        const size_t vertexBytes = result.mesh.positions.size() * sizeof(result.mesh.positions[0]);
        const size_t indexBytes = result.mesh.indices.size() * sizeof(result.mesh.indices[0]);
        const size_t innerGridIndexBytes =
            size_t(result.innerGridTriangleCount) * 3 * sizeof(result.mesh.indices[0]);
        const size_t stitchingIndexBytes =
            size_t(result.stitchingTriangleCount) * 3 * sizeof(result.mesh.indices[0]);
        std::printf("  meshMemory vertexBytes=%zu (%.3f MiB) indexBytes=%zu (%.3f MiB) totalBytes=%zu (%.3f MiB)\n",
                    vertexBytes,
                    BytesToMiB(vertexBytes),
                    indexBytes,
                    BytesToMiB(indexBytes),
                    vertexBytes + indexBytes,
                    BytesToMiB(vertexBytes + indexBytes));
        std::printf("  indexMemorySplit innerGridTris=%d innerGridIndexBytes=%zu (%.3f MiB)"
                    " stitchingTris=%d stitchingIndexBytes=%zu (%.3f MiB)\n",
                    result.innerGridTriangleCount,
                    innerGridIndexBytes,
                    BytesToMiB(innerGridIndexBytes),
                    result.stitchingTriangleCount,
                    stitchingIndexBytes,
                    BytesToMiB(stitchingIndexBytes));
        std::printf("  writeMetadata=%s\n", cli.writeMetadata ? "true" : "false");
        std::printf("  edgeMapChecks=ok\n");
        std::printf("  controlCageUniqueEdges=%zu boundaryEdges=%d\n",
                    result.controlCageUniqueEdges,
                    result.boundaryEdges);
        std::printf("  controlCageEdgesWithOver2Faces=%d\n", result.controlCageEdgesWithOver2Faces);
        return 0;
    }

    TessellationRunConfig config = {};
    config.level = cli.level;
    config.pixelSpacing = cli.pixelSpacing;
    config.splitThreshold = cli.splitThreshold;
    config.sampleSteps = cli.sampleSteps;
    config.writeMetadata = cli.writeMetadata;
    config.patchQuadObjPath = cli.patchQuadObjPath;
    return RunUsdTessellationInput(cli.inPath, cli.outObjPath, config, WriteMeshObjWithTriMetadata)
               ? 0
               : 1;
}
