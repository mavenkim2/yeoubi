#pragma once

#include "io/usd/load.h"
#include "scene/scene.h"
#include "tessellation/subdivision.h"
#include "util/assert.h"
#include "util/float3.h"

#include <pxr/usd/usd/primRange.h>
#include <pxr/usd/usd/stage.h>
#include <pxr/usd/usdGeom/camera.h>
#include <pxr/usd/usdGeom/mesh.h>
#include <pxr/usd/usdGeom/xformCache.h>

#include <cctype>
#include <cstdio>
#include <filesystem>
#include <string>
#include <vector>

struct TessellationRunConfig
{
    int level = 3;
    float pixelSpacing = 1.0f;
    int splitThreshold = 1;
    int sampleSteps = 8;
    bool writeMetadata = false;
    std::string patchQuadObjPath;
};

struct TessellationOutputEntry
{
    size_t index = 0;
    std::string primPath;
    std::string outObjPath;
    size_t verts = 0;
    size_t tris = 0;
};

using TessellationMeshWriterFn = bool (*)(const ybi::Mesh &mesh,
                                          const std::vector<int> &triPatchFaceIds,
                                          const std::vector<int> &triCoarseFaceIds,
                                          const std::vector<int> &triPtexFaceIds,
                                          const std::vector<int> &triQuadrants,
                                          const std::string &outObjPath,
                                          bool writeMetadata);

static std::string Lowercase(const std::string &s)
{
    std::string out = s;
    for (char &c : out)
    {
        c = char(std::tolower(static_cast<unsigned char>(c)));
    }
    return out;
}

static bool IsJsonInputPath(const std::string &path)
{
    const std::string lower = Lowercase(path);
    return lower.size() >= 5 && lower.substr(lower.size() - 5) == ".json";
}

static bool IsUsdInputPath(const std::string &path)
{
    const std::string lower = Lowercase(path);
    return (lower.size() >= 4 && lower.substr(lower.size() - 4) == ".usd") ||
           (lower.size() >= 5 && lower.substr(lower.size() - 5) == ".usda") ||
           (lower.size() >= 5 && lower.substr(lower.size() - 5) == ".usdc");
}

static ybi::float3 ComputeSubdivisionMeshCenter(const ybi::SubdivisionMesh &mesh)
{
    if (mesh.vertices.size() == 0)
    {
        return ybi::make_float3(0.0f);
    }
    ybi::float3 sum = ybi::make_float3(0.0f);
    for (const ybi::float3 &p : mesh.vertices)
    {
        sum += p;
    }
    return sum * (1.0f / float(mesh.vertices.size()));
}

static bool GetUsdCameraAndSubdivPrimPaths(const std::string &usdPath,
                                           ybi::float3 *outCameraWorldPos,
                                           std::string *outCameraPath,
                                           std::vector<std::string> *outSubdivPrimPaths)
{
    pxr::UsdStageRefPtr stage = pxr::UsdStage::Open(usdPath);
    if (!stage)
    {
        return false;
    }

    if (outSubdivPrimPaths)
    {
        outSubdivPrimPaths->clear();
    }

    pxr::UsdGeomXformCache xformCache(pxr::UsdTimeCode::Default());
    const pxr::Usd_PrimFlagsConjunction flags =
        pxr::UsdPrimIsActive && pxr::UsdPrimIsLoaded && !pxr::UsdPrimIsAbstract;
    const pxr::Usd_PrimFlagsPredicate predicate(flags);

    bool foundCamera = false;
    for (const pxr::UsdPrim &prim : stage->Traverse(predicate))
    {
        if (!foundCamera && prim.IsA<pxr::UsdGeomCamera>())
        {
            const pxr::GfMatrix4d localToWorld = xformCache.GetLocalToWorldTransform(prim);
            const pxr::GfVec3d p = localToWorld.Transform(pxr::GfVec3d(0.0, 0.0, 0.0));
            if (outCameraWorldPos)
            {
                *outCameraWorldPos = ybi::make_float3(float(p[0]), float(p[1]), float(p[2]));
            }
            if (outCameraPath)
            {
                *outCameraPath = prim.GetPath().GetString();
            }
            foundCamera = true;
        }

        if (outSubdivPrimPaths && prim.IsA<pxr::UsdGeomMesh>())
        {
            pxr::UsdGeomMesh mesh(prim);
            pxr::TfToken scheme;
            if (mesh.GetSubdivisionSchemeAttr().Get(&scheme) && scheme == pxr::UsdGeomTokens->catmullClark)
            {
                outSubdivPrimPaths->push_back(prim.GetPath().GetString());
            }
        }
    }

    return foundCamera;
}

static std::string SanitizePrimPathForFilename(const std::string &primPath)
{
    std::string out;
    out.reserve(primPath.size() + 8);
    for (char c : primPath)
    {
        if ((c >= 'a' && c <= 'z') || (c >= 'A' && c <= 'Z') || (c >= '0' && c <= '9') || c == '_')
        {
            out.push_back(c);
        }
        else
        {
            out.push_back('_');
        }
    }
    return out.empty() ? "mesh" : out;
}

static std::string JsonEscape(const std::string &s)
{
    std::string out;
    out.reserve(s.size() + 8);
    for (char c : s)
    {
        if (c == '\\' || c == '"')
        {
            out.push_back('\\');
        }
        out.push_back(c);
    }
    return out;
}

static bool WriteTessellationManifest(const std::filesystem::path &manifestPath,
                                      const std::vector<TessellationOutputEntry> &entries)
{
    FILE *f = std::fopen(manifestPath.string().c_str(), "w");
    if (!f)
    {
        return false;
    }

    std::fprintf(f, "{\n  \"outputs\": [\n");
    for (size_t i = 0; i < entries.size(); ++i)
    {
        const TessellationOutputEntry &e = entries[i];
        std::fprintf(f,
                     "    {\"meshIndex\": %zu, \"primPath\": \"%s\", \"obj\": \"%s\", \"verts\": %zu, \"tris\": %zu}%s\n",
                     e.index,
                     JsonEscape(e.primPath).c_str(),
                     JsonEscape(e.outObjPath).c_str(),
                     e.verts,
                     e.tris,
                     (i + 1 < entries.size()) ? "," : "");
    }
    std::fprintf(f, "  ]\n}\n");
    std::fclose(f);
    return true;
}

static bool RunUsdTessellationInput(const std::string &inPath,
                                    const std::string &outPathString,
                                    const TessellationRunConfig &config,
                                    TessellationMeshWriterFn meshWriter)
{
    if (!meshWriter)
    {
        std::fprintf(stderr, "Internal error: null mesh writer\n");
        return false;
    }

    ybi::ScenePool scenePool;
    ybi::LoadUSDScene(&scenePool, inPath);
    if (scenePool.scenes.empty() || scenePool.rootSceneIndex >= scenePool.scenes.size())
    {
        std::fprintf(stderr, "Failed to load USD scene: %s\n", inPath.c_str());
        return false;
    }

    ybi::float3 cameraWorldPos = ybi::make_float3(0.0f);
    std::string cameraPath;
    std::vector<std::string> subdivPrimPaths;
    const bool hasCamera =
        GetUsdCameraAndSubdivPrimPaths(inPath, &cameraWorldPos, &cameraPath, &subdivPrimPaths);
    YBI_ASSERT(hasCamera);
    if (!hasCamera)
    {
        std::fprintf(stderr, "USD camera missing: %s\n", inPath.c_str());
        return false;
    }

    struct MeshRef
    {
        const ybi::SubdivisionMesh *mesh = nullptr;
        std::string primPath;
    };
    std::vector<MeshRef> meshRefs;
    size_t primPathCursor = 0;
    for (const std::unique_ptr<ybi::Scene> &scenePtr : scenePool.scenes)
    {
        if (!scenePtr)
        {
            continue;
        }
        for (const ybi::SubdivisionMesh &mesh : scenePtr->subdivisionMeshes)
        {
            const std::string primPath =
                primPathCursor < subdivPrimPaths.size() ? subdivPrimPaths[primPathCursor] : "";
            meshRefs.push_back({&mesh, primPath});
            primPathCursor++;
        }
    }

    if (meshRefs.empty())
    {
        std::fprintf(stderr, "No subdivision meshes found in USD: %s\n", inPath.c_str());
        return false;
    }

    const bool singleOutput = meshRefs.size() == 1;
    const std::filesystem::path outPath(outPathString);
    if (singleOutput)
    {
        const std::filesystem::path parent = outPath.parent_path();
        if (!parent.empty())
        {
            std::filesystem::create_directories(parent);
        }
    }
    else
    {
        if (std::filesystem::exists(outPath) && !std::filesystem::is_directory(outPath))
        {
            std::fprintf(stderr, "Output must be a directory for multiple meshes: %s\n", outPathString.c_str());
            return false;
        }
        std::filesystem::create_directories(outPath);
    }

    std::filesystem::path patchQuadDir;
    if (!config.patchQuadObjPath.empty() && !singleOutput)
    {
        patchQuadDir = std::filesystem::path(config.patchQuadObjPath);
        if (std::filesystem::exists(patchQuadDir) && !std::filesystem::is_directory(patchQuadDir))
        {
            std::fprintf(stderr,
                         "patchquadobj must be a directory for multiple meshes: %s\n",
                         config.patchQuadObjPath.c_str());
            return false;
        }
        std::filesystem::create_directories(patchQuadDir);
    }

    std::vector<TessellationOutputEntry> outputs;
    outputs.reserve(meshRefs.size());
    size_t totalVerts = 0;
    size_t totalTris = 0;

    for (size_t i = 0; i < meshRefs.size(); ++i)
    {
        const MeshRef &ref = meshRefs[i];
        YBI_ASSERT(ref.mesh);

        std::string objPath = outPathString;
        std::string patchQuadPath;
        if (!singleOutput)
        {
            char prefix[32];
            std::snprintf(prefix, sizeof(prefix), "%04zu", i);
            const std::string safePrim = SanitizePrimPathForFilename(ref.primPath);
            objPath = (outPath / (std::string(prefix) + "_" + safePrim + ".obj")).string();
            if (!config.patchQuadObjPath.empty())
            {
                patchQuadPath =
                    (patchQuadDir / (std::string(prefix) + "_" + safePrim + "_patch_quads.obj")).string();
            }
        }
        else if (!config.patchQuadObjPath.empty())
        {
            patchQuadPath = config.patchQuadObjPath;
        }

        ybi::SubdivisionRunOptions options = {};
        options.level = config.level;
        options.pixelSpacing = config.pixelSpacing;
        options.splitThreshold = config.splitThreshold;
        options.sampleSteps = config.sampleSteps;
        options.eye = cameraWorldPos;
        options.lookAt = ComputeSubdivisionMeshCenter(*ref.mesh);
        options.subdivisionScheme = "catmullClark";
        options.patchQuadObjPath = patchQuadPath;

        ybi::SubdivisionRunResult result = {};
        if (!ybi::SubdivideAdaptive(*ref.mesh, options, &result))
        {
            std::fprintf(stderr, "Subdivision failed for mesh %zu\n", i);
            return false;
        }

        if (!meshWriter(result.mesh,
                        result.trianglePatchFaceIds,
                        result.triangleCoarseFaceIds,
                        result.trianglePtexFaceIds,
                        result.triangleQuadrants,
                        objPath,
                        config.writeMetadata))
        {
            std::fprintf(stderr, "Failed to write OBJ: %s\n", objPath.c_str());
            return false;
        }

        outputs.push_back({i,
                           ref.primPath.empty() ? std::string("<unknown>") : ref.primPath,
                           objPath,
                           result.mesh.positions.size(),
                           result.mesh.indices.size() / 3});
        totalVerts += result.mesh.positions.size();
        totalTris += result.mesh.indices.size() / 3;

        std::printf("  mesh[%zu] prim=%s out=%s verts=%zu tris=%zu\n",
                    i,
                    outputs.back().primPath.c_str(),
                    objPath.c_str(),
                    outputs.back().verts,
                    outputs.back().tris);
    }

    if (!singleOutput)
    {
        const std::filesystem::path manifestPath = outPath / "tessellation_manifest.json";
        if (!WriteTessellationManifest(manifestPath, outputs))
        {
            std::fprintf(stderr, "Failed to write manifest: %s\n", manifestPath.string().c_str());
            return false;
        }
        std::printf("  manifest=%s\n", manifestPath.string().c_str());
    }
    std::printf("  camera=%s eye=(%.6f %.6f %.6f)\n",
                cameraPath.c_str(),
                double(cameraWorldPos.x),
                double(cameraWorldPos.y),
                double(cameraWorldPos.z));
    std::printf("  meshes=%zu totalVerts=%zu totalTris=%zu writeMetadata=%s\n",
                outputs.size(),
                totalVerts,
                totalTris,
                config.writeMetadata ? "true" : "false");
    return true;
}
