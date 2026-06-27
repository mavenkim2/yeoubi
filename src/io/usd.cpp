#include "io/usd.h"
#include "core/affine_transform.h"
#include "core/allocation.h"
#include "core/arena.h"
#include "core/array.h"
#include "core/assert.h"
#include "io/usd.h"
#include "scene/scene.h"
#include <cmath>
#include <filesystem>
#include <memory>
#include <pxr/base/gf/camera.h>
#include <pxr/base/gf/frustum.h>
#include <pxr/base/vt/array.h>
#include <pxr/usd/usd/prim.h>
#include <pxr/usd/usd/primFlags.h>
#include <pxr/usd/usd/stage.h>
#include <pxr/usd/usdGeom/basisCurves.h>
#include <pxr/usd/usdGeom/camera.h>
#include <pxr/usd/usdGeom/mesh.h>
#include <pxr/usd/usdGeom/pointInstancer.h>
#include <pxr/usd/usdGeom/tokens.h>
#include <pxr/usd/usdGeom/xformCache.h>
#include <pxr/usd/usdShade/material.h>
#include <pxr/usd/usdShade/materialBindingAPI.h>
#include <string>
#include <tbb/parallel_for.h>
#include <unordered_set>
#include <vector>

namespace ybi
{

#define USD_YBI_ASSERT(expr) expr

struct USDPrimTraversalData
{
    pxr::UsdPrim prim;
    std::vector<pxr::UsdPrim> meshes;
    std::vector<pxr::UsdPrim> pointInstancers;
    std::vector<pxr::UsdPrim> instances;
    std::vector<pxr::UsdPrim> curves;

    std::vector<pxr::UsdGeomCamera> cameras;
};

static pxr::UsdShadeMaterial GetPrimMaterial(const pxr::UsdPrim &prim)
{
    pxr::UsdShadeMaterialBindingAPI bindingApi(prim);
    pxr::UsdShadeMaterial material = bindingApi.ComputeBoundMaterial(pxr::UsdShadeTokens->full);
    if (!material)
    {
        material = bindingApi.ComputeBoundMaterial(pxr::UsdShadeTokens->preview);
    }
    if (!material)
    {
        material = bindingApi.ComputeBoundMaterial(pxr::UsdShadeTokens->allPurpose);
    }

    if (material)
    {
        // printf("Prim %s is bound to material: %s\n",
        //        prim.GetPath().GetText(),
        //        material.GetPath().GetText());
    }
    else
    {
        // printf("Prim %s has no material bound.\n", prim.GetPath().GetText());
    }

    return material;
}

std::unique_ptr<SceneLoader> SceneLoader::CreateSceneLoader(const std::string &filePath)
{
    std::string ext = std::filesystem::path(filePath).extension().string();
    if (ext.find("usd") != std::string::npos)
    {
        std::unique_ptr<SceneLoader> loader = std::make_unique<USDSceneLoader>(filePath);
        return loader;
    }
    else
    {
        return nullptr;
    }
}

USDSceneLoader::USDSceneLoader(std::string filePath) : filePath(filePath)
{
    allowedPurposes.push_back(pxr::UsdGeomTokens->default_);
    allowedPurposes.push_back(pxr::UsdGeomTokens->render);
}

AffineTransform USDSceneLoader::ConvertAffineTransform(const pxr::GfMatrix4d &m) const
{
    const pxr::GfMatrix4d t = m.GetTranspose();
    const pxr::GfVec4d r0 = t.GetRow(0);
    const pxr::GfVec4d r1 = t.GetRow(1);
    const pxr::GfVec4d r2 = t.GetRow(2);

    return AffineTransform((float)r0[0],
                           (float)r0[1],
                           (float)r0[2],
                           (float)r0[3],
                           (float)r1[0],
                           (float)r1[1],
                           (float)r1[2],
                           (float)r1[3],
                           (float)r2[0],
                           (float)r2[1],
                           (float)r2[2],
                           (float)r2[3]);
}

Matrix4x4 USDSceneLoader::ConvertMatrix4x4(const pxr::GfMatrix4d &m) const
{
    const pxr::GfMatrix4d t = m.GetTranspose();
    const pxr::GfVec4d r0 = t.GetRow(0);
    const pxr::GfVec4d r1 = t.GetRow(1);
    const pxr::GfVec4d r2 = t.GetRow(2);
    const pxr::GfVec4d r3 = t.GetRow(3);

    return Matrix4x4((float)r0[0],
                     (float)r0[1],
                     (float)r0[2],
                     (float)r0[3],
                     (float)r1[0],
                     (float)r1[1],
                     (float)r1[2],
                     (float)r1[3],
                     (float)r2[0],
                     (float)r2[1],
                     (float)r2[2],
                     (float)r2[3],
                     (float)r3[0],
                     (float)r3[1],
                     (float)r3[2],
                     (float)r3[3]);
}

void USDSceneLoader::ConvertMeshToTriangles(const pxr::VtIntArray &faceIndices,
                                            const pxr::VtIntArray &faceCounts,
                                            Array<int> &triangles) const
{
    size_t index = 0;
    for (int faceCount : faceCounts)
    {
        YBI_ASSERT(faceCount >= 3);

        for (size_t i = 1; i < faceCount - 1; i++)
        {
            triangles.EmplaceBack(faceIndices[index]);
            triangles.EmplaceBack(faceIndices[index + i]);
            triangles.EmplaceBack(faceIndices[index + i + 1]);
        }

        index += faceCount;
    }
}

bool USDSceneLoader::IsPurposeAllowed(pxr::UsdPrim prim) const
{
    pxr::UsdGeomImageable imageable(prim);

    pxr::TfToken purpose = imageable.ComputePurpose();
    for (const pxr::TfToken &allowed : allowedPurposes)
    {
        if (purpose == allowed)
        {
            return true;
        }
    }
    return false;
}

void USDSceneLoader::ProcessPolygonMesh(Scene *scene,
                                        pxr::UsdPrim sceneRoot,
                                        size_t meshOffset,
                                        size_t meshIndex,
                                        const pxr::UsdPrim prim,
                                        pxr::UsdGeomXformCache &cache,
                                        pxr::UsdTimeCode timeCode) const
{
    pxr::UsdGeomMesh usdMesh(prim);

    pxr::UsdGeomImageable imageable(prim);
    // prim.ComputePurpo

    pxr::VtVec3fArray usdPositions;
    pxr::VtIntArray faceIndices;
    pxr::VtIntArray faceCounts;
    pxr::TfToken subdivisionScheme;

    bool result = usdMesh.GetPointsAttr().Get(&usdPositions, timeCode);
    YBI_ASSERT(result);
    result = usdMesh.GetFaceVertexIndicesAttr().Get(&faceIndices, timeCode);
    YBI_ASSERT(result);
    result = usdMesh.GetFaceVertexCountsAttr().Get(&faceCounts, timeCode);
    YBI_ASSERT(result);
    result = usdMesh.GetSubdivisionSchemeAttr().Get(&subdivisionScheme, timeCode);
    YBI_ASSERT(result);

    std::string topologyError;
    if (!pxr::UsdGeomMesh::ValidateTopology(
            faceIndices, faceCounts, usdPositions.size(), &topologyError))
    {
        fprintf(stderr,
                "Invalid USD mesh topology: %s\nPrim: %s\n",
                topologyError.c_str(),
                prim.GetPath().GetText());
        YBI_ASSERT(false);
    }

    Array<Vec3f> positions(usdPositions.size());
    memcpy(positions.data(), usdPositions.data(), sizeof(Vec3f) * usdPositions.size());
    Array<int> triangleIndices;
    ConvertMeshToTriangles(faceIndices, faceCounts, triangleIndices);

    Mesh mesh(std::move(positions), std::move(triangleIndices));

    size_t index = meshOffset + meshIndex;
    scene->container->meshes[index] = std::move(mesh);

    bool resetsXformStack = false;
    pxr::GfMatrix4d localFromMesh =
        cache.ComputeRelativeTransform(prim, sceneRoot, &resetsXformStack);
    Instance &instance = scene->instances[meshIndex];
    instance.transform = ConvertAffineTransform(localFromMesh);
    instance.handle.type = GeometryType::Mesh;
    instance.handle.index = index;
}

void USDSceneLoader::TraversePrim(USDPrimTraversalData *data,
                                  pxr::UsdPrim root,
                                  bool skipPrototypePrims) const
{
    const pxr::Usd_PrimFlagsPredicate filterPredicate =
        pxr::UsdPrimIsActive && pxr::UsdPrimIsLoaded && !pxr::UsdPrimIsAbstract;

    std::vector<pxr::UsdPrim> stack;
    stack.push_back(root);

    while (!stack.empty())
    {
        const pxr::UsdPrim prim = stack.back();
        stack.pop_back();

        if (skipPrototypePrims && (prim.IsInPrototype() || prim.IsInstanceProxy()))
        {
            continue;
        }

        pxr::UsdGeomImageable imageable(prim);
        if (imageable && imageable.ComputeVisibility() == pxr::UsdGeomTokens->invisible)
        {
            continue;
        }

        if (!IsPurposeAllowed(prim))
        {
            continue;
        }

        bool pushChildren = true;
        if (prim.IsInstance())
        {
            data->instances.push_back(prim);
            pushChildren = false;
        }
        else if (prim.IsA<pxr::UsdGeomPointInstancer>())
        {
            data->pointInstancers.push_back(prim);
            pushChildren = false;
        }
        else if (prim.IsA<pxr::UsdGeomMesh>())
        {
            data->meshes.push_back(prim);
        }
        else if (prim.IsA<pxr::UsdGeomBasisCurves>())
        {
            data->curves.push_back(prim);
        }
        else if (prim.IsA<pxr::UsdGeomCamera>())
        {
            data->cameras.push_back(pxr::UsdGeomCamera(prim));
        }

        if (!pushChildren)
        {
            continue;
        }

        for (const pxr::UsdPrim &child : prim.GetFilteredChildren(filterPredicate))
        {
            stack.push_back(child);
        }
    }
}

static void ProcessCatmullClarkMesh(DataContainer *container,
                                    pxr::UsdPrim meshPrim,
                                    Mesh *mesh,
                                    pxr::UsdTimeCode timeCode = pxr::UsdTimeCode::Default())
{
    pxr::VtVec3fArray positions;
    pxr::VtIntArray faceIndices;
    pxr::VtIntArray faceCounts;
    pxr::VtVec3fArray normals;

    pxr::TfToken interpolationBoundary;
    pxr::TfToken fvLinearInterpolation;
    pxr::TfToken triangleSubdivisionRule;

    pxr::VtIntArray cornerIndices;
    pxr::VtFloatArray cornerSharpnesses;
    pxr::VtIntArray creaseIndices;
    pxr::VtIntArray creaseLengths;
    pxr::VtFloatArray creaseSharpnesses;
    pxr::VtIntArray holeIndices;

    pxr::UsdGeomMesh usdMesh(meshPrim);

    USD_YBI_ASSERT(usdMesh.GetPointsAttr().Get(&positions, timeCode));

    USD_YBI_ASSERT(usdMesh.GetFaceVertexIndicesAttr().Get(&faceIndices, timeCode));
    USD_YBI_ASSERT(usdMesh.GetFaceVertexCountsAttr().Get(&faceCounts, timeCode));

    if (usdMesh.GetNormalsAttr().Get(&normals, timeCode))
    {
        pxr::TfToken normalsInterpToken = usdMesh.GetNormalsInterpolation();
        for (auto n : normals)
        {
            // printf("help: %f %f %f\n", n[0], n[1], n[2]);
        }
    }

    USD_YBI_ASSERT(
        usdMesh.GetTriangleSubdivisionRuleAttr().Get(&triangleSubdivisionRule, timeCode));
    USD_YBI_ASSERT(usdMesh.GetInterpolateBoundaryAttr().Get(&interpolationBoundary, timeCode));
    USD_YBI_ASSERT(
        usdMesh.GetFaceVaryingLinearInterpolationAttr().Get(&fvLinearInterpolation, timeCode));

    usdMesh.GetCornerIndicesAttr().Get(&cornerIndices, timeCode);
    usdMesh.GetCornerSharpnessesAttr().Get(&cornerSharpnesses, timeCode);
    usdMesh.GetCreaseIndicesAttr().Get(&creaseIndices, timeCode);
    usdMesh.GetCreaseLengthsAttr().Get(&creaseLengths, timeCode);
    usdMesh.GetCreaseSharpnessesAttr().Get(&creaseSharpnesses, timeCode);
    usdMesh.GetHoleIndicesAttr().Get(&holeIndices, timeCode);

    pxr::UsdPrim prim = usdMesh.GetPrim();

    // mesh->SetNumVertices(positions.size());

#if 0
    Array<Vec3> positionsArray(positions);
    Array<int> faceIndicesArray(faceIndices);
    Array<int> faceCountsArray(faceCounts);

    BoundaryInterpolation interpolation = BOUNDARY_INTERPOLATION_EDGE_AND_CORNER;
    if (interpolationBoundary == "none")
    {
        interpolation = BOUNDARY_INTERPOLATION_NONE;
    }
    else if (interpolationBoundary == "edgeOnly")
    {
        interpolation = BOUNDARY_INTERPOLATION_EDGE;
    }
    else if (interpolationBoundary == "edgeAndCorner")
    {
        interpolation = BOUNDARY_INTERPOLATION_EDGE_AND_CORNER;
    }

    FVarLinearInterpolation fvarLinear = FVAR_LINEAR_CORNERS_ONLY;
    if (fvLinearInterpolation == "none")
    {
        fvarLinear = FVAR_LINEAR_NONE;
    }
    else if (fvLinearInterpolation == "cornersOnly")
    {
        fvarLinear = FVAR_LINEAR_CORNERS_ONLY;
    }
    else if (fvLinearInterpolation == "cornersPlus1")
    {
        fvarLinear = FVAR_LINEAR_CORNERS_PLUS1;
    }
    else if (fvLinearInterpolation == "cornersPlus2")
    {
        fvarLinear = FVAR_LINEAR_CORNERS_PLUS2;
    }
    else if (fvLinearInterpolation == "boundaries")
    {
        fvarLinear = FVAR_LINEAR_BOUNDARIES;
    }
    else if (fvLinearInterpolation == "all")
    {
        fvarLinear = FVAR_LINEAR_ALL;
    }

    TriangleSubdivisionRule triangleSubdivision = TRIANGLE_SUBDIVISION_CATMULL_CLARK;
    if (triangleSubdivisionRule == "smooth")
    {
        triangleSubdivision = TRIANGLE_SUBDIVISION_SMOOTH;
    }

    Array<int> cornerIndicesArray(cornerIndices);
    Array<float> cornerSharpnessesArray(cornerSharpnesses);
    Array<int> creaseIndicesArray(creaseIndices);
    Array<int> creaseLengthsArray(creaseLengths);
    Array<float> creaseSharpnessesArray(creaseSharpnesses);

    scene->subdivisionMeshes.emplace_back(std::move(positionsArray),
                                          std::move(faceIndicesArray),
                                          std::move(faceCountsArray),
                                          std::move(cornerIndicesArray),
                                          std::move(cornerSharpnessesArray),
                                          std::move(creaseIndicesArray),
                                          std::move(creaseLengthsArray),
                                          std::move(creaseSharpnessesArray),
                                          std::move(holeIndices),
                                          interpolation,
                                          fvarLinear,
                                          triangleSubdivision);
    SubdivisionMesh &subdivMesh = scene->subdivisionMeshes.back();
    subdivMesh.primPath = prim.GetPath().GetString();
    subdivMesh.materialIndex = materialIndex;
    subdivMesh.parentFromLocal = parentFromLocal;
    Attribute stAttr;
    if (BuildSubdivisionMeshSTAttribute(mesh, faceIndices, &stAttr, timeCode))
    {
        subdivMesh.attributes.push_back(std::move(stAttr));
    }
#endif
}

static void PrintMatrix(const pxr::GfMatrix4d &d)
{
    printf("%f %f %f %f\n % f % f % f % f\n % f % f % f % f\n % f % f % f % f\n ",
           d[0][0],
           d[0][1],
           d[0][2],
           d[0][3],
           d[1][0],
           d[1][1],
           d[1][2],
           d[1][3],
           d[2][0],
           d[2][1],
           d[2][2],
           d[2][3],
           d[3][0],
           d[3][1],
           d[3][2],
           d[3][3]);
}

static bool RoughlyEqual(double a, double b)
{
    return std::abs(b - a) < 1e-6;
}

static bool RoughlyEqual(const pxr::GfMatrix4d &a, const pxr::GfMatrix4d &b)
{
    for (size_t row = 0; row < 4; row++)
    {
        for (size_t column = 0; column < 4; column++)
        {
            if (!RoughlyEqual(a[row][column], b[row][column]))
            {
                return false;
            }
        }
    }

    return true;
}

Scene *USDSceneLoader::LoadScene(Arena *arena, DataContainer *container)
{
    pxr::UsdStageRefPtr stage = pxr::UsdStage::Open(filePath.c_str());
    pxr::UsdPrim root = stage->GetPseudoRoot();

    USDPrimTraversalData collectiveData;

    // std::vector<pxr::UsdPrim> traversablePrims;
    std::unordered_map<std::string, size_t> prototypeToIndexMap;
    size_t prototypeCount = 1;

    std::vector<USDPrimTraversalData> traversalDatas;
    traversalDatas.push_back({root});
    bool skipPrototypePrims = true;

    size_t traversablePrimIndex = 0;

    while (traversablePrimIndex != traversalDatas.size())
    {
        size_t index = traversablePrimIndex++;
        USDPrimTraversalData &data = traversalDatas[index];
        TraversePrim(&data, data.prim, skipPrototypePrims);

        skipPrototypePrims = false;

        size_t numInstances = data.instances.size();
        size_t numPointInstancers = data.pointInstancers.size();

        for (size_t instanceIndex = 0; instanceIndex < numInstances; instanceIndex++)
        {
            pxr::UsdPrim instance = traversalDatas[index].instances[instanceIndex];
            const pxr::UsdPrim prototype = instance.GetPrototype();
            YBI_ASSERT(prototype);

            const std::string path = prototype.GetPath().GetString();

            // Prune duplicates
            if (prototypeToIndexMap.find(path) == prototypeToIndexMap.end())
            {
                prototypeToIndexMap.insert({path, prototypeCount++});
                traversalDatas.push_back({prototype});
                YBI_ASSERT(prototypeCount == traversalDatas.size());
            }
        }

        for (size_t instanceIndex = 0; instanceIndex < numPointInstancers; instanceIndex++)
        {
            pxr::UsdPrim instance = traversalDatas[index].pointInstancers[instanceIndex];
            pxr::UsdGeomPointInstancer pointInstancer(instance);
            pxr::SdfPathVector prototypePathsForInstancer;
            bool result =
                pointInstancer.GetPrototypesRel().GetTargets(&prototypePathsForInstancer);

            YBI_ASSERT(result);

            for (const pxr::SdfPath &path : prototypePathsForInstancer)
            {
                std::string pathString = path.GetString();
                if (prototypeToIndexMap.find(pathString) == prototypeToIndexMap.end())
                {
                    prototypeToIndexMap.insert({pathString, prototypeCount++});
                    pxr::UsdPrim proto = stage->GetPrimAtPath(path);
                    traversalDatas.push_back({proto});
                    YBI_ASSERT(prototypeCount == traversalDatas.size());
                }
            }
        }
    }

    printf("num scenes: %lli\n", traversalDatas.size());

    const pxr::UsdTimeCode timeCode(stage->GetStartTimeCode());
    pxr::UsdGeomXformCache xformCache(timeCode);

    // Load cameras
    USDPrimTraversalData &data = traversalDatas[0];
    YBI_ASSERT(!data.cameras.empty());

    printf("num cameras: %zu\n", data.cameras.size());
    pxr::UsdGeomCamera camera = data.cameras[0];
    pxr::UsdPrim cameraPrim = camera.GetPrim();
    pxr::GfCamera gfCam = camera.GetCamera(timeCode);
    pxr::GfFrustum frustum = gfCam.GetFrustum();
    pxr::GfMatrix4d viewMatrix = frustum.ComputeViewMatrix();
    pxr::GfMatrix4d projMatrix = frustum.ComputeProjectionMatrix();

    container->camera.cameraFromWorld = ConvertAffineTransform(viewMatrix);
    container->camera.clipFromCamera = ConvertMatrix4x4(projMatrix);

    // TODO: cli option
    const int targetWidth = 1920;
    const int targetHeight = 1080;

    float sx = targetWidth * 0.5f;
    float sy = targetHeight * 0.5f;
    Matrix4x4 rasterFromNDC = Matrix4x4(sx, 0, 0, sx, 0, -sy, 0, sy, 0, 0, 1, 0, 0, 0, 0, 1);

    container->camera.rasterFromCamera = rasterFromNDC * container->camera.clipFromCamera;
    container->camera.width = targetWidth;
    container->camera.height = targetHeight;

    std::vector<Scene> scenes;
    scenes.resize(traversalDatas.size());

    std::vector<pxr::UsdShadeMaterial> materials;

    tbb::enumerable_thread_specific<pxr::UsdGeomXformCache> xformCaches(
        [&]() { return pxr::UsdGeomXformCache(timeCode); });

    for (size_t i = 0; i < traversalDatas.size(); i++)
    {
        Scene *scene = &scenes[i];
        scene->container = container;

        USDPrimTraversalData &data = traversalDatas[i];
        size_t totalInstanceCount = data.instances.size();
        for (pxr::UsdPrim &prim : data.pointInstancers)
        {
            pxr::UsdGeomPointInstancer pointInstancer(prim);
            totalInstanceCount += pointInstancer.GetInstanceCount(timeCode);
        }

        size_t meshOffset = scene->container->meshes.size();
        size_t dstInstanceIndex = data.meshes.size();
        scene->instances.Resize(totalInstanceCount + data.meshes.size());
        scene->container->meshes.resize(meshOffset + data.meshes.size());

        tbb::parallel_for(tbb::blocked_range<size_t>(0, data.meshes.size()),
                          [&](const tbb::blocked_range<size_t> &range) {
                              pxr::UsdGeomXformCache &cache = xformCaches.local();
                              for (size_t i = range.begin(); i != range.end(); i++)
                              {
                                  pxr::UsdPrim prim = data.meshes[i];

                                  ProcessPolygonMesh(
                                      scene, data.prim, meshOffset, i, prim, cache, timeCode);
                              }
                          });

        // Get all materials
        for (size_t i = 0; i < data.meshes.size(); i++)
        {
            pxr::UsdPrim prim = data.meshes[i];
            pxr::UsdShadeMaterial material = GetPrimMaterial(prim);

            if (!material)
            {
                printf("%s is missing a material\n", prim.GetPath().GetString().c_str());
                continue;
            }
            materials.push_back(material);
        }

        for (pxr::UsdPrim &usdInstance : data.instances)
        {
            pxr::UsdPrim prototype = usdInstance.GetPrototype();

            auto it = prototypeToIndexMap.find(prototype.GetPath().GetString());
            YBI_ASSERT(it != prototypeToIndexMap.end());

            bool resetsXformStack = false;
            const pxr::GfMatrix4d instanceTransform =
                xformCache.ComputeRelativeTransform(usdInstance, data.prim, &resetsXformStack);

            Instance &instance = scene->instances[dstInstanceIndex++];
            instance.transform = ConvertAffineTransform(instanceTransform);
            instance.handle.type = GeometryType::Scene;
            instance.handle.index = it->second;

            YBI_ASSERT(!(instance.transform == AffineTransform::Zero()));
            YBI_ASSERT(it->second != 0);
        }

        for (pxr::UsdPrim &prim : data.pointInstancers)
        {
            pxr::UsdGeomPointInstancer pointInstancer(prim);

            pxr::VtIntArray protoIndices;
            pxr::SdfPathVector prototypePaths;

            bool result = pointInstancer.GetProtoIndicesAttr().Get(&protoIndices, timeCode);
            YBI_ASSERT(result);
            result = pointInstancer.GetPrototypesRel().GetTargets(&prototypePaths);
            YBI_ASSERT(result);
            std::vector<bool> mask = pointInstancer.ComputeMaskAtTime(timeCode);

            bool resetsXformStack = false;
            const pxr::GfMatrix4d pointInstancerLocalToScene =
                xformCache.ComputeRelativeTransform(prim, data.prim, &resetsXformStack);

            pxr::VtArray<pxr::GfMatrix4d> instanceTransforms;
            result = pointInstancer.ComputeInstanceTransformsAtTime(
                &instanceTransforms,
                timeCode,
                timeCode,
                pxr::UsdGeomPointInstancer::IncludeProtoXform,
                pxr::UsdGeomPointInstancer::ApplyMask);
            YBI_ASSERT(result);

            result = pxr::UsdGeomPointInstancer::ApplyMaskToArray(mask, &protoIndices);
            YBI_ASSERT(result);

            YBI_ASSERT(instanceTransforms.size() == protoIndices.size());

            for (size_t i = 0; i < instanceTransforms.size(); i++)
            {
                int protoIndex = protoIndices[i];
                YBI_ASSERT(protoIndex >= 0);
                YBI_ASSERT(static_cast<size_t>(protoIndex) < prototypePaths.size());

                const pxr::SdfPath &prototypePath = prototypePaths[protoIndex];
                auto it = prototypeToIndexMap.find(prototypePath.GetString());
                YBI_ERROR(it != prototypeToIndexMap.end(),
                          "help: %s\n",
                          prototypePath.GetString().c_str());

                const pxr::GfMatrix4d parentFromLocal =
                    instanceTransforms[i] * pointInstancerLocalToScene;

                YBI_ASSERT(!(ConvertAffineTransform(parentFromLocal) == AffineTransform::Zero()));

                Instance &instance = scene->instances[dstInstanceIndex++];
                instance.transform = ConvertAffineTransform(parentFromLocal);
                instance.handle.type = GeometryType::Scene;
                instance.handle.index = it->second;

                YBI_ASSERT(it->second != 0);
            }
        }

        scene->instances.Resize(dstInstanceIndex);
    }

    // Handle materials
    for (pxr::UsdShadeMaterial &material : materials)
    {
        // MaterialInfo info = {};
        // info.materialPath = material.GetPath().GetString();
        // info.packed = {};

        if (pxr::UsdShadeShader shader = material.ComputeSurfaceSource())
        {
            pxr::TfToken token;
            if (shader.GetShaderId(&token) && token == pxr::TfToken("UsdPreviewSurface"))
            {
                for (const pxr::UsdShadeInput &input : shader.GetInputs())
                {
                    // std::string texturePath;
                    // std::string swizzle;
                    // TextureWrapMode wrapS = TEXTURE_WRAP_MODE_UNKNOWN;
                    // TextureWrapMode wrapT = TEXTURE_WRAP_MODE_UNKNOWN;

                    auto sources = input.GetConnectedSources();

                    if (sources.size() == 0)
                        continue;
                    YBI_ERROR(sources.size() == 1, "Num sources: %zu\n", sources.size());

                    pxr::UsdPrim sourcePrim = sources[0].source.GetPrim();
                    YBI_ASSERT(sourcePrim.IsA<pxr::UsdShadeShader>());
                    pxr::UsdShadeShader shader(sourcePrim);

                    pxr::TfToken shaderId;
                    shader.GetShaderId(&shaderId);
                    YBI_ASSERT(shaderId == pxr::TfToken("UsdUVTexture"));

                    pxr::UsdShadeInput fileInput = shader.GetInput(pxr::TfToken("file"));
                    YBI_ASSERT(fileInput);

                    // if (!TryGetImageTexturePath(input, texturePath, &swizzle, &wrapS, &wrapT))
                    // {
                    //     continue;
                    // }

                    // MaterialTextureInput textureInput = {};
                    // textureInput.inputName = input.GetBaseName().GetString();
                    // textureInput.texturePath = texturePath;
                    // textureInput.swizzle = swizzle;
                    // textureInput.wrapS = wrapS;
                    // textureInput.wrapT = wrapT;
                    // outInfo->textureInputs.push_back(std::move(textureInput));
                }
            }
            else
            {
                YBI_ASSERT(0);
                // TODO: support non-UsdPreviewSurface material networks.
                // FinalizePackedMaterial(&info.packed);
            }
        }
        else if (pxr::UsdShadeShader shader = material.ComputeSurfaceSource(pxr::TfToken("mdl")))
        {
            YBI_ASSERT(0);
        }
        else
        {
            // NOTE: empty materials???
            for (const pxr::UsdShadeOutput &output : material.GetOutputs())
            {
                printf("material output: %s\n", output.GetFullName().GetText());
            }
            // YBI_LOGFATAL("%s\n", material.GetPath().GetString().c_str());
            // YBI_ASSERT(0);
        }
    }

    // Flatten
    Scene *newScene = new Scene(container);

    std::vector<std::pair<Scene *, AffineTransform>> stack;
    stack.push_back(std::make_pair(&scenes[0], AffineTransform::Identity()));

    while (!stack.empty())
    {
        std::pair<Scene *, AffineTransform> info = stack.back();
        stack.pop_back();

        size_t index = 0;
        for (Instance &instance : info.first->instances)
        {
            AffineTransform transform = info.second * instance.transform;

            switch (instance.handle.type)
            {
                case GeometryType::Scene:
                {
                    stack.push_back(std::make_pair(&scenes[instance.handle.index], transform));
                    break;
                }
                case GeometryType::Mesh:
                {
                    newScene->instances.EmplaceBack(transform, instance.handle);
                    break;
                }
                default:
                {
                    // TODO: not yet implemented
                    YBI_ASSERT(0);
                    break;
                }
            }
            index++;
        }
    }

    printf("num instances: %lli\n", newScene->instances.size());

    // NOTE: by this point all instances should be traversed
#if 0
    container->meshes = arena->Alloc<Mesh>(data.meshes.size());
    for (size_t meshIndex = 0; meshIndex < data.meshes.size(); meshIndex++)
    {
        ProcessCatmullClarkMesh(container, data.meshes[meshIndex], &container->meshes[meshIndex]);
        printf("num vertices: %lli\n", container->meshes[meshIndex].GetNumVertices());
    }
#endif

    return newScene;
}

} // namespace ybi
