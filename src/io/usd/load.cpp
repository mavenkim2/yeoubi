#include "io/usd/load.h"
#include "io/usd/instance_dag_build.h"
#include "pxr/base/gf/vec2f.h"
#include "scene/attributes.h"
#include "scene/scene.h"
#include "util/assert.h"
#include "util/float2.h"
#include "util/float3.h"
#include "util/float3x4.h"
#include "util/float4.h"
#include "util/float4x4.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <pxr/base/gf/frustum.h>
#include <pxr/base/gf/matrix4d.h>
#include <pxr/base/gf/matrix4f.h>
#include <pxr/base/gf/vec4f.h>
#include <pxr/base/vt/types.h>
#include <pxr/base/vt/value.h>
#include <pxr/usd/sdf/path.h>
#include <pxr/usd/sdf/assetPath.h>
#include <pxr/usd/sdf/valueTypeName.h>
#include <pxr/usd/usd/prim.h>
#include <pxr/usd/usd/primRange.h>
#include <pxr/usd/usd/relationship.h>
#include <pxr/usd/usd/stage.h>
#include <pxr/usd/usd/variantSets.h>
#include <pxr/usd/usdGeom/basisCurves.h>
#include <pxr/usd/usdGeom/camera.h>
#include <pxr/usd/usdGeom/curves.h>
#include <pxr/usd/usdGeom/hermiteCurves.h>
#include <pxr/usd/usdGeom/mesh.h>
#include <pxr/usd/usdGeom/metrics.h>
#include <pxr/usd/usdGeom/nurbsCurves.h>
#include <pxr/usd/usdGeom/pointInstancer.h>
#include <pxr/usd/usdGeom/primvarsAPI.h>
#include <pxr/usd/usdGeom/scope.h>
#include <pxr/usd/usdGeom/subset.h>
#include <pxr/usd/usdGeom/xformCache.h>
#include <pxr/usd/usdGeom/xform.h>
#include <pxr/usd/usdGeom/xformable.h>
#include <pxr/usd/usdRender/settings.h>
#include <pxr/usd/usdShade/material.h>
#include <pxr/usd/usdShade/materialBindingAPI.h>
#include <pxr/usd/usdShade/shader.h>
#include <pxr/usd/usdVol/volume.h>
#include <unordered_map>
#include <utility>

#include <immintrin.h>

YBI_NAMESPACE_BEGIN

#define USD_ASSERT(expr) \
    { \
        bool result = expr; \
        if (!result) \
        { \
            printf("USD Error in %s (%s:%d)\n", #expr, __FILE__, __LINE__); \
            YBI_ASSERT(false); \
        } \
    }

static void CollectUSDCameras(const pxr::UsdPrim &root, std::vector<pxr::UsdGeomCamera> *out)
{
    YBI_ASSERT(out);
    out->clear();

    const pxr::Usd_PrimFlagsPredicate filterPredicate =
        pxr::UsdPrimIsActive && pxr::UsdPrimIsLoaded && !pxr::UsdPrimIsAbstract;
    std::vector<pxr::UsdPrim> stack;
    stack.push_back(root);
    while (!stack.empty())
    {
        const pxr::UsdPrim prim = stack.back();
        stack.pop_back();

        if (prim.IsA<pxr::UsdGeomCamera>())
        {
            out->push_back(pxr::UsdGeomCamera(prim));
        }

        for (const pxr::UsdPrim &child : prim.GetFilteredChildren(filterPredicate))
        {
            stack.push_back(child);
        }
    }
}

static bool
CreateRuntimeScenesFromBuildSceneDAG(const USDBuildSceneDAG &dag,
                                     ScenePool *scenePool,
                                     std::string *error)
{
    YBI_ASSERT(scenePool);
    if (dag.rootSceneIndex >= dag.scenes.size())
    {
        if (error)
        {
            *error = "build scene dag root index out of range";
        }
        return false;
    }

    scenePool->scenes.clear();
    scenePool->scenes.reserve(dag.scenes.size());
    for (size_t i = 0; i < dag.scenes.size(); i++)
    {
        scenePool->scenes.push_back(std::make_unique<Scene>());
    }

    for (uint32_t sceneIndex = 0; sceneIndex < dag.scenes.size(); sceneIndex++)
    {
        const USDBuildScene &buildScene = dag.scenes[sceneIndex];
        Scene *outScene = scenePool->scenes[sceneIndex].get();
        outScene->meshes.clear();
        outScene->curves.clear();
        outScene->instances.clear();
        outScene->childScenes.clear();
        outScene->meshes.reserve(buildScene.meshes.size());
        outScene->curves.reserve(buildScene.curves.size());
        outScene->instances.reserve(buildScene.instances.size());
        outScene->childScenes.reserve(buildScene.instances.size());

        for (const USDBuildSceneInstance &instance : buildScene.instances)
        {
            if (instance.childSceneIndex >= dag.scenes.size())
            {
                if (error)
                {
                    *error = "child scene index out of range while converting build scenes";
                }
                return false;
            }
            outScene->childScenes.push_back(scenePool->scenes[instance.childSceneIndex].get());
            outScene->instances.emplace_back(
                instance.parentFromLocal,
                static_cast<uint32_t>(outScene->childScenes.size() - 1));
        }
    }

    scenePool->rootSceneIndex = dag.rootSceneIndex;
    return true;
}

static pxr::UsdShadeMaterial GetPrimMaterial(const pxr::UsdPrim &prim,
                                             pxr::TfToken token = pxr::UsdShadeTokens->full)
{
    pxr::UsdShadeMaterialBindingAPI bindingApi(prim);
    // TODO: set this as an option
    pxr::UsdShadeMaterial material = bindingApi.ComputeBoundMaterial(token);

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

static int AddMaterialToMap(std::unordered_map<std::string, int> &materialMap,
                            std::vector<pxr::UsdShadeMaterial> &materials,
                            const pxr::UsdShadeMaterial &material)
{
    int materialIndex = -1;
    if (material)
    {
        pxr::SdfPath sdfPath = material.GetPath();
        if (!sdfPath.IsEmpty())
        {
            std::string path = sdfPath.GetString();

            auto found = materialMap.find(path);
            if (found == materialMap.end())
            {
                materialIndex = static_cast<int>(materials.size());
                materialMap.emplace(path, materialIndex);
                materials.push_back(material);
            }
            else
            {
                materialIndex = found->second;
            }
        }
    }
    return materialIndex;
}

// NOTE: For shade inputs that most likely have a single UsdUVTexture
static bool TryGetSingleConnectedShader(const pxr::UsdShadeInput &input,
                                        pxr::UsdShadeShader &shaderOut)
{
    auto sources = input.GetConnectedSources();
    if (sources.size() != 1)
    {
        // TODO: handle multiple or zero sources for a shading input.
        printf("Texture gather: unsupported source count (%zu) for input %s\n",
               sources.size(),
               input.GetBaseName().GetText());
        return false;
    }

    pxr::UsdPrim sourcePrim = sources[0].source.GetPrim();
    if (!sourcePrim.IsA<pxr::UsdShadeShader>())
    {
        // TODO: handle non-shader source nodes (e.g. nodegraphs/other source types).
        printf("Texture gather: unsupported source prim type %s for input %s\n",
               sourcePrim.GetTypeName().GetText(),
               input.GetBaseName().GetText());
        return false;
    }
    shaderOut = pxr::UsdShadeShader(sourcePrim);
    return true;
}

static std::string OutputNameToSwizzle(const pxr::TfToken &sourceName)
{
    const std::string name = sourceName.GetString();
    if (name == "r")
    {
        return "R";
    }
    if (name == "g")
    {
        return "G";
    }
    if (name == "b")
    {
        return "B";
    }
    if (name == "a")
    {
        return "A";
    }
    if (name == "rgb")
    {
        return "RGB";
    }
    if (name == "rgba")
    {
        return "RGBA";
    }
    return "";
}

static bool TryGetImageTexturePath(const pxr::UsdShadeInput &input,
                                   std::string &outPath,
                                   std::string *outSwizzle)
{
    outPath.clear();
    if (outSwizzle)
    {
        outSwizzle->clear();
    }

    pxr::UsdShadeShader sourceShader;
    if (!TryGetSingleConnectedShader(input, sourceShader))
    {
        return false;
    }

    pxr::TfToken shaderId;
    sourceShader.GetShaderId(&shaderId);
    if (shaderId != pxr::TfToken("UsdUVTexture"))
    {
        // TODO: support non-UsdUVTexture image producers (MaterialX image nodes, etc.).
        printf("Texture gather: unsupported shader node type %s at %s\n",
               shaderId.GetText(),
               sourceShader.GetPath().GetText());
        return false;
    }

    pxr::UsdShadeInput fileInput = sourceShader.GetInput(pxr::TfToken("file"));
    if (!fileInput)
    {
        // TODO: support alternate file-bearing texture nodes/ports.
        printf("Texture gather: UsdUVTexture missing 'file' input at %s\n",
               sourceShader.GetPath().GetText());
        return false;
    }

    pxr::SdfAssetPath assetPath;
    if (!fileInput.Get(&assetPath))
    {
        // TODO: support connected/indirected file inputs.
        printf("Texture gather: could not read 'file' value at %s\n",
               sourceShader.GetPath().GetText());
        return false;
    }

    outPath = assetPath.GetResolvedPath();
    if (outPath.empty())
    {
        outPath = assetPath.GetAssetPath();
    }
    if (outPath.empty())
    {
        printf("Texture gather: empty texture file path at %s\n", sourceShader.GetPath().GetText());
        return false;
    }

    if (outSwizzle)
    {
        auto sources = input.GetConnectedSources();
        if (sources.size() == 1)
        {
            outSwizzle->assign(OutputNameToSwizzle(sources[0].sourceName));
        }
    }

    return true;
}

static void ReadNtcMaterialInfo(const pxr::UsdShadeMaterial &material,
                                ScenePool::MaterialInfo *outInfo)
{
    YBI_ASSERT(outInfo);
    outInfo->ntcDiffuseFile.clear();
    outInfo->ntcDiffuseTextureName.clear();

    const pxr::UsdPrim materialPrim = material.GetPrim();
    if (!materialPrim)
    {
        return;
    }

    const pxr::UsdAttribute ntcFileAttr =
        materialPrim.GetAttribute(pxr::TfToken("ybi:ntc:diffuseFile"));
    if (ntcFileAttr)
    {
        pxr::SdfAssetPath assetPath;
        if (ntcFileAttr.Get(&assetPath))
        {
            outInfo->ntcDiffuseFile = assetPath.GetResolvedPath();
            if (outInfo->ntcDiffuseFile.empty())
            {
                outInfo->ntcDiffuseFile = assetPath.GetAssetPath();
            }
        }
        else
        {
            // Back-compat if authored as plain string.
            ntcFileAttr.Get(&outInfo->ntcDiffuseFile);
        }
    }

    const pxr::UsdAttribute ntcTextureNameAttr =
        materialPrim.GetAttribute(pxr::TfToken("ybi:ntc:diffuseTextureName"));
    if (ntcTextureNameAttr)
    {
        ntcTextureNameAttr.Get(&outInfo->ntcDiffuseTextureName);
    }
    if (!outInfo->ntcDiffuseFile.empty() && outInfo->ntcDiffuseTextureName.empty())
    {
        outInfo->ntcDiffuseTextureName = "diffuseColor";
    }
}

static void AppendUniqueTexturePath(std::vector<std::string> &paths, const std::string &path)
{
    if (std::find(paths.begin(), paths.end(), path) == paths.end())
    {
        paths.push_back(path);
    }
}

static int GetMaterialIndex(const std::unordered_map<std::string, int> &materialMap,
                            const pxr::UsdShadeMaterial &material)
{
    if (!material)
    {
        return -1;
    }
    const std::string materialPath = material.GetPath().GetString();
    auto it = materialMap.find(materialPath);
    if (it == materialMap.end())
    {
        return -1;
    }
    return it->second;
}

static bool BuildTriangulatedMeshTexcoords(const pxr::UsdGeomMesh &mesh,
                                           const pxr::VtIntArray &faceCounts,
                                           const pxr::VtIntArray &faceIndices,
                                           int numTriangles,
                                           Array<float2> *outTexcoords,
                                           Array<int> *outTexcoordIndices)
{
    YBI_ASSERT(outTexcoords);
    YBI_ASSERT(outTexcoordIndices);

    pxr::UsdGeomPrimvarsAPI primvars(mesh);
    pxr::UsdGeomPrimvar stPrimvar = primvars.GetPrimvar(pxr::TfToken("st"));
    if (!stPrimvar)
    {
        return false;
    }

    pxr::VtVec2fArray stValues;
    if (!stPrimvar.Get(&stValues, 0.0))
    {
        return false;
    }
    if (stValues.empty())
    {
        return false;
    }

    pxr::VtIntArray stIndices;
    const bool hasExplicitIndices = stPrimvar.IsIndexed() && stPrimvar.GetIndices(&stIndices, 0.0);
    const pxr::TfToken interpolation = stPrimvar.GetInterpolation();

    Array<float2> texcoords(stValues);
    Array<int> triTexcoordIndices(static_cast<size_t>(numTriangles) * 3u);

    int faceIndexOffset = 0;
    int triIndexOffset = 0;
    for (int faceCount : faceCounts)
    {
        int faceTc[4] = {-1, -1, -1, -1};
        if (faceCount < 3 || faceCount > 4)
        {
            return false;
        }

        for (int corner = 0; corner < faceCount; corner++)
        {
            const int cornerIndex = faceIndexOffset + corner;
            int tcIndex = -1;
            if (hasExplicitIndices)
            {
                if (cornerIndex < 0 || cornerIndex >= static_cast<int>(stIndices.size()))
                {
                    return false;
                }
                tcIndex = stIndices[cornerIndex];
            }
            else if (interpolation == pxr::UsdGeomTokens->faceVarying)
            {
                tcIndex = cornerIndex;
            }
            else
            {
                if (cornerIndex < 0 || cornerIndex >= static_cast<int>(faceIndices.size()))
                {
                    return false;
                }
                tcIndex = faceIndices[cornerIndex];
            }

            if (tcIndex < 0 || tcIndex >= static_cast<int>(texcoords.size()))
            {
                return false;
            }
            faceTc[corner] = tcIndex;
        }

        if (faceCount == 3)
        {
            triTexcoordIndices[triIndexOffset++] = faceTc[0];
            triTexcoordIndices[triIndexOffset++] = faceTc[1];
            triTexcoordIndices[triIndexOffset++] = faceTc[2];
        }
        else
        {
            triTexcoordIndices[triIndexOffset++] = faceTc[0];
            triTexcoordIndices[triIndexOffset++] = faceTc[1];
            triTexcoordIndices[triIndexOffset++] = faceTc[2];
            triTexcoordIndices[triIndexOffset++] = faceTc[0];
            triTexcoordIndices[triIndexOffset++] = faceTc[2];
            triTexcoordIndices[triIndexOffset++] = faceTc[3];
        }
        faceIndexOffset += faceCount;
    }

    outTexcoords->Resize(texcoords.size());
    if (texcoords.size() > 0)
    {
        memcpy(outTexcoords->data(), texcoords.data(), sizeof(float2) * texcoords.size());
    }
    outTexcoordIndices->Resize(triTexcoordIndices.size());
    if (triTexcoordIndices.size() > 0)
    {
        memcpy(outTexcoordIndices->data(),
               triTexcoordIndices.data(),
               sizeof(int) * triTexcoordIndices.size());
    }
    return true;
}

static bool BuildSubdivisionMeshTexcoords(const pxr::UsdGeomMesh &mesh,
                                          const pxr::VtIntArray &faceIndices,
                                          Array<float2> *outTexcoords,
                                          Array<int> *outTexcoordIndices)
{
    YBI_ASSERT(outTexcoords);
    YBI_ASSERT(outTexcoordIndices);

    pxr::UsdGeomPrimvarsAPI primvars(mesh);
    pxr::UsdGeomPrimvar stPrimvar = primvars.GetPrimvar(pxr::TfToken("st"));
    if (!stPrimvar)
    {
        return false;
    }

    pxr::VtVec2fArray stValues;
    if (!stPrimvar.Get(&stValues, 0.0))
    {
        return false;
    }
    if (stValues.empty())
    {
        return false;
    }

    pxr::VtIntArray stIndices;
    const bool hasExplicitIndices = stPrimvar.IsIndexed() && stPrimvar.GetIndices(&stIndices, 0.0);
    const pxr::TfToken interpolation = stPrimvar.GetInterpolation();

    Array<float2> texcoords(stValues);
    Array<int> cornerTexcoordIndices(faceIndices.size());
    for (size_t corner = 0; corner < faceIndices.size(); ++corner)
    {
        int tcIndex = -1;
        if (hasExplicitIndices)
        {
            if (corner >= stIndices.size())
            {
                return false;
            }
            tcIndex = stIndices[corner];
        }
        else if (interpolation == pxr::UsdGeomTokens->faceVarying)
        {
            tcIndex = int(corner);
        }
        else if (interpolation == pxr::UsdGeomTokens->vertex ||
                 interpolation == pxr::UsdGeomTokens->varying)
        {
            tcIndex = faceIndices[corner];
        }
        else
        {
            return false;
        }

        if (tcIndex < 0 || tcIndex >= int(texcoords.size()))
        {
            return false;
        }
        cornerTexcoordIndices[corner] = tcIndex;
    }

    *outTexcoords = std::move(texcoords);
    *outTexcoordIndices = std::move(cornerTexcoordIndices);
    return true;
}

static PrimvarInterpolation ConvertPrimvarInterpolation(pxr::TfToken &token)
{
    if (token == "constant")
    {
        return PrimvarInterpolation::Constant;
    }
    if (token == "uniform")
    {
        return PrimvarInterpolation::Uniform;
    }
    if (token == "vertex")
    {
        return PrimvarInterpolation::Vertex;
    }
    if (token == "varying")
    {
        return PrimvarInterpolation::Varying;
    }
    if (token == "faceVarying")
    {
        return PrimvarInterpolation::FaceVarying;
    }
    return PrimvarInterpolation::Unknown;
}

static AttributeType ConvertPrimvarTypeName(pxr::SdfValueTypeName type)
{
    if (type == pxr::SdfValueTypeNames->FloatArray)
    {
        return AttributeType::Float;
    }
    if (type == pxr::SdfValueTypeNames->Double)
    {
        return AttributeType::Float;
    }
    if (type == pxr::SdfValueTypeNames->UCharArray)
    {
        return AttributeType::Int8;
    }
    if (type == pxr::SdfValueTypeNames->IntArray || type == pxr::SdfValueTypeNames->Int)
    {
        return AttributeType::Int;
    }
    if (type == pxr::SdfValueTypeNames->Float2Array)
    {
        return AttributeType::Float2;
    }
    if (type == pxr::SdfValueTypeNames->TexCoord2dArray)
    {
        return AttributeType::Float2;
    }
    if (type == pxr::SdfValueTypeNames->TexCoord2fArray)
    {
        return AttributeType::Float2;
    }
    if (type == pxr::SdfValueTypeNames->TexCoord2hArray)
    {
        return AttributeType::Float2;
    }
    if (type == pxr::SdfValueTypeNames->TexCoord3dArray)
    {
        return AttributeType::Float3;
    }
    if (type == pxr::SdfValueTypeNames->TexCoord3fArray)
    {
        return AttributeType::Float3;
    }
    if (type == pxr::SdfValueTypeNames->TexCoord3hArray)
    {
        return AttributeType::Float3;
    }
    if (type == pxr::SdfValueTypeNames->Float3Array)
    {
        return AttributeType::Float3;
    }
    if (type == pxr::SdfValueTypeNames->Point3fArray)
    {
        return AttributeType::Float3;
    }
    if (type == pxr::SdfValueTypeNames->Point3dArray)
    {
        return AttributeType::Float3;
    }
    if (type == pxr::SdfValueTypeNames->Point3hArray)
    {
        return AttributeType::Float3;
    }
    if (type == pxr::SdfValueTypeNames->Normal3fArray)
    {
        return AttributeType::Float3;
    }
    if (type == pxr::SdfValueTypeNames->Normal3dArray)
    {
        return AttributeType::Float3;
    }
    if (type == pxr::SdfValueTypeNames->Normal3hArray)
    {
        return AttributeType::Float3;
    }
    if (type == pxr::SdfValueTypeNames->Vector3fArray)
    {
        return AttributeType::Float3;
    }
    if (type == pxr::SdfValueTypeNames->Vector3hArray)
    {
        return AttributeType::Float3;
    }
    if (type == pxr::SdfValueTypeNames->Vector3dArray)
    {
        return AttributeType::Float3;
    }
    if (type == pxr::SdfValueTypeNames->Color3fArray)
    {
        return AttributeType::Float3;
    }
    if (type == pxr::SdfValueTypeNames->Color3hArray)
    {
        return AttributeType::Float3;
    }
    if (type == pxr::SdfValueTypeNames->Color3dArray)
    {
        return AttributeType::Float3;
    }
    if (type == pxr::SdfValueTypeNames->Color4fArray)
    {
        return AttributeType::Float4;
    }
    if (type == pxr::SdfValueTypeNames->Color4hArray)
    {
        return AttributeType::Float4;
    }
    if (type == pxr::SdfValueTypeNames->Color4dArray)
    {
        return AttributeType::Float4;
    }
    if (type == pxr::SdfValueTypeNames->BoolArray || type == pxr::SdfValueTypeNames->Bool)
    {
        return AttributeType::Bool;
    }
    if (type == pxr::SdfValueTypeNames->QuatfArray)
    {
        return AttributeType::Quaternion;
    }
    if (type == pxr::SdfValueTypeNames->QuatdArray)
    {
        return AttributeType::Quaternion;
    }
    if (type == pxr::SdfValueTypeNames->QuathArray)
    {
        return AttributeType::Quaternion;
    }
    if (type == pxr::SdfValueTypeNames->Token)
    {
        return AttributeType::Unknown;
    }
    YBI_LOGFATAL("Unsupported primvar type: %s", type.GetAsToken().GetText());
    return AttributeType::Unknown;
}

static MemoryView<uint8_t> ConvertPrimvarValues(Scene *scene, const pxr::VtValue &values)
{
    if (values.IsHolding<pxr::VtArray<pxr::GfVec2f>>())
    {
        pxr::VtArray<pxr::GfVec2f> array = values.Get<pxr::VtArray<pxr::GfVec2f>>();
        // MemoryView<float2> view = scene->arena.PushArray<float2>(array.size());
        // memcpy(view.data(), array.data(), array.size() * sizeof(float2));
        // return view.CastToBytes();
    }
    else
    {
        YBI_LOGFATAL("Unsupported primvar values type: %s", values.GetTypeName().c_str());
    }

    return MemoryView<uint8_t>();
}

static void ProcessPrimvars(pxr::UsdPrim prim, pxr::UsdTimeCode timeCode, Scene *scene)
{
    pxr::UsdGeomPrimvarsAPI primvarsAPI(prim);
    for (const pxr::UsdGeomPrimvar &primVar : primvarsAPI.GetPrimvars())
    {
        pxr::VtValue values;
        if (primVar.HasValue())
        {
            bool flattened = primVar.IsIndexed();
            USD_ASSERT(primVar.ComputeFlattened(&values, timeCode));
            pxr::TfToken interpolationToken = primVar.GetInterpolation();
            PrimvarInterpolation interpolation = ConvertPrimvarInterpolation(interpolationToken);

            pxr::SdfValueTypeName typeName = primVar.GetTypeName();
            AttributeType attrType = ConvertPrimvarTypeName(typeName);

            pxr::TfToken name = pxr::UsdGeomPrimvar::StripPrimvarsName(primVar.GetName());

            if (name == "st")
            {
                MemoryView<uint8_t> data = ConvertPrimvarValues(scene, values);
                scene->attributes.EmplaceBack(data, attrType, interpolation);
            }
            // printf("name: %s, typeName: %s, interp: %s\n",
            //        name.GetText(),
            //        typeName.GetAsToken().GetText(),
            //        interpolationToken.GetText());
        }
    }
}

static void
ProcessCatmullClarkMesh(pxr::UsdGeomMesh &mesh, Scene *scene, pxr::UsdTimeCode timeCode = 0.0)
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

    USD_ASSERT(mesh.GetPointsAttr().Get(&positions, timeCode));

    USD_ASSERT(mesh.GetFaceVertexIndicesAttr().Get(&faceIndices, timeCode));
    USD_ASSERT(mesh.GetFaceVertexCountsAttr().Get(&faceCounts, timeCode));

    if (mesh.GetNormalsAttr().Get(&normals, timeCode))
    {
        pxr::TfToken normalsInterpToken = mesh.GetNormalsInterpolation();
        for (auto n : normals)
        {
            // printf("help: %f %f %f\n", n[0], n[1], n[2]);
        }
    }

    USD_ASSERT(mesh.GetTriangleSubdivisionRuleAttr().Get(&triangleSubdivisionRule, timeCode));
    USD_ASSERT(mesh.GetInterpolateBoundaryAttr().Get(&interpolationBoundary, timeCode));
    USD_ASSERT(mesh.GetFaceVaryingLinearInterpolationAttr().Get(&fvLinearInterpolation, timeCode));

    mesh.GetCornerIndicesAttr().Get(&cornerIndices, timeCode);
    mesh.GetCornerSharpnessesAttr().Get(&cornerSharpnesses, timeCode);
    mesh.GetCreaseIndicesAttr().Get(&creaseIndices, timeCode);
    mesh.GetCreaseLengthsAttr().Get(&creaseLengths, timeCode);
    mesh.GetCreaseSharpnessesAttr().Get(&creaseSharpnesses, timeCode);
    mesh.GetHoleIndicesAttr().Get(&holeIndices, timeCode);

    pxr::UsdPrim prim = mesh.GetPrim();
    size_t attributeStart = scene->attributes.size();
    ProcessPrimvars(prim, timeCode, scene);
    size_t attributeEnd = scene->attributes.size();

    Array<float3> positionsArray(positions);
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
                                          attributeStart,
                                          attributeEnd,
                                          interpolation,
                                          fvarLinear,
                                          triangleSubdivision);
    SubdivisionMesh &subdivMesh = scene->subdivisionMeshes.back();
    subdivMesh.primPath = prim.GetPath().GetString();
    BuildSubdivisionMeshTexcoords(mesh, faceIndices, &subdivMesh.texcoords, &subdivMesh.texcoordIndices);
}

static void ProcessUSDBasisCurve(pxr::UsdGeomBasisCurves &curve, Scene *scene)
{
    size_t numCurves = curve.GetCurveCount(0.0);

    pxr::VtIntArray curveVertexCounts;
    pxr::VtVec3fArray points;
    pxr::VtFloatArray widths;
    pxr::TfToken basisToken;
    pxr::TfToken typeToken;
    pxr::TfToken wrapToken;

    USD_ASSERT(curve.GetCurveVertexCountsAttr().Get(&curveVertexCounts, 0.0));
    USD_ASSERT(curve.GetPointsAttr().Get(&points, 0.0));
    USD_ASSERT(curve.GetWidthsAttr().Get(&widths, 0.0));
    USD_ASSERT(curve.GetBasisAttr().Get(&basisToken, 0.0));
    USD_ASSERT(curve.GetTypeAttr().Get(&typeToken, 0.0));
    USD_ASSERT(curve.GetWrapAttr().Get(&wrapToken, 0.0));

    uint32_t offset = 0;
    uint32_t curveIndex = 0;

#if 0
    uint32_t numBits = 0;
    for (auto count : curveVertexCounts)
    {
        numBits += 32 * 3;

        float minF[3];
        float maxF[3];
        uint32_t expMin[3];
        uint32_t expMax[3];
        uint32_t mantissaMin[3];
        uint32_t mantissaMax[3];
        for (int i = 0; i < 3; i++)
        {
            minF[i] = std::numeric_limits<float>::infinity();
            maxF[i] = -std::numeric_limits<float>::infinity();

            expMin[i] = 512;
            expMax[i] = 0u;

            mantissaMin[i] = ~0u;
            mantissaMax[i] = 0;
        }
        for (uint32_t i = offset; i < offset + count; i++)
        {
            auto point = points[i];

            for (int j = 0; j < 3; j++)
            {
                minF[j] = std::min(minF[j], point[j]);
                maxF[j] = std::max(maxF[j], point[j]);

                uint32_t test;
                memcpy(&test, &point[j], sizeof(float));

                uint32_t exponent = (test >> 23) & 0xff;
                uint32_t mantissa = test & 0x7fffff;
                if (exponent == 0)
                {
                    // printf("0 exp: %f %f %f, %u\n", point[0], point[1], point[2], test);
                }
                expMin[j] = std::min(expMin[j], exponent);
                expMax[j] = std::max(expMax[j], exponent);

                mantissaMin[j] = std::min(mantissaMin[j], mantissa);
                mantissaMax[j] = std::max(mantissaMax[j], mantissa);
            }
        }

        if (numCurves == 3689911)
        {
            printf("curve: %u\n", curveIndex++);
            for (uint32_t i = offset; i < offset + count; i++)
            {
                auto point = points[i];
                printf("point: %f %f %f\n", point[0], point[1], point[2]);
            }
        }

        numBits += 8 * 3 + 23 * 3;
        for (int j = 0; j < 3; j++)
        {
            float delta = maxF[j] - minF[j];
            int expDelta = expMax[j] - expMin[j];
            int mantissaDelta = mantissaMax[j] - mantissaMin[j];

            int numExpBits, numMantissaBits;
            if (expDelta == 0)
            {
                numExpBits = 0;
            }
            else
            {
                _BitScanReverse((unsigned long *)&numExpBits, expDelta);
                numExpBits++;
            }

            if (mantissaDelta == 0)
            {
                numMantissaBits = 0;
            }
            else
            {
                _BitScanReverse((unsigned long *)&numMantissaBits, mantissaDelta);
                numMantissaBits++;
            }
            numBits += count * (numExpBits + numMantissaBits);

            int test;
            memcpy(&test, &delta, sizeof(float));

            // printf("exponents: %i %i %i\n", expMin[j], expMax[j], expDelta);
            // printf("mant: %i %i %i\n", mantissaMin[j], mantissaMax[j], mantissaDelta);
            // printf("test delta: floats %u %f %f\n", test, maxF[j], delta);
        }
        offset += count;
    }
    printf("uncompressed: %i, compressed: %i\n", points.size() * 12, numBits / 8);
#endif

    std::vector<pxr::UsdGeomSubset> subsets = pxr::UsdGeomSubset::GetAllGeomSubsets(curve);
    if (subsets.size())
    {
        printf("has subsets\n");
    }

    printf("curve: %s %s %s\n", basisToken.GetText(), typeToken.GetText(), wrapToken.GetText());

    int vStep = 1;

    if (basisToken.GetString() == "bezier")
    {
        vStep = 3;
    }

    int totalNumVertices = 0;

    Array<float3> positions(points);
    Array<float> curveWidths(widths);
    Array<int> curveOffsets(numCurves);
    int curveOffsetIndex = 0;

    for (int curveVertexCount : curveVertexCounts)
    {
        curveOffsets[curveOffsetIndex++] = totalNumVertices;
        totalNumVertices += curveVertexCount;
    }

    YBI_ASSERT(totalNumVertices == points.size());
    scene->curves.emplace_back(
        std::move(positions), std::move(curveWidths), std::move(curveOffsets));

    int curveFlags = 0;

    if (typeToken == "cubic")
    {
        curveFlags |= CurveFlags::CURVE_FLAGS_CUBIC;
    }
    else
    {
        YBI_ASSERT(typeToken == "linear");
        curveFlags |= CurveFlags::CURVE_FLAGS_LINEAR;
    }

    if (curve.GetNormalsAttr().HasValue())
    {
        pxr::VtVec3fArray normals;

        USD_ASSERT(curve.GetNormalsAttr().Get(&normals, 0.0));
        printf("normal %f %f %f\n", normals[0][0], normals[0][1], normals[0][2]);
        curveFlags |= CurveFlags::CURVE_FLAGS_RIBBON;
    }
    else
    {
        curveFlags |= CurveFlags::CURVE_FLAGS_TUBE;
    }

    printf("basis: %zi %zi %s\n", numCurves, points.size(), basisToken.GetText());
}

void LoadUSDScene(ScenePool *scenePool, const std::string &filePath)
{
    YBI_ASSERT(scenePool);
    scenePool->camera = Camera();
    pxr::UsdStageRefPtr stage = pxr::UsdStage::Open(filePath.c_str());

    if (!stage)
    {
        printf("error opening usd stage at %s\n", filePath.c_str());
        return;
    }

    if (pxr::UsdGeomGetStageUpAxis(stage) != pxr::UsdGeomTokens->z)
    {
    }

    double startTimeCode = stage->GetStartTimeCode();
    double endTimeCode = stage->GetEndTimeCode();
    double fps = stage->GetFramesPerSecond();
    double tcps = stage->GetTimeCodesPerSecond();
    double timeCodesPerFrame = tcps / fps;

    printf("start: %f, end: %f, fps: %f, tcps: %f\n", startTimeCode, endTimeCode, fps, tcps);

    pxr::UsdRenderSettings settings = pxr::UsdRenderSettings::GetStageRenderSettings(stage);
    if (settings)
    {
        printf("found\n");
        float shutterOpen, shutterClose;
        if (settings.GetPrim().GetAttribute(pxr::TfToken("shutter:open")).Get(&shutterOpen))
        {
            printf("Global Render Shutter Open: %f\n", shutterOpen);
        }
    }

    pxr::UsdPrim root = stage->GetPseudoRoot();
    USDBuildSceneDAG buildSceneDAG = {};
    std::string buildSceneError;
    if (!BuildInstanceDAGFromUSD(stage, &buildSceneDAG, &buildSceneError))
    {
        printf("build instance dag failed: %s\n", buildSceneError.c_str());
        return;
    }
    std::vector<pxr::UsdGeomCamera> cameras;
    CollectUSDCameras(root, &cameras);

    std::string createScenesError;
    if (!CreateRuntimeScenesFromBuildSceneDAG(buildSceneDAG, scenePool, &createScenesError))
    {
        printf("build scene creation failed: %s\n", createScenesError.c_str());
        return;
    }

    if (cameras.size())
    {
        pxr::UsdGeomCamera &camera = cameras[0];
        const pxr::UsdTimeCode timeCode(0.0);
        pxr::GfCamera gfCam = camera.GetCamera(timeCode);
        pxr::GfFrustum frustum = gfCam.GetFrustum();
        pxr::GfMatrix4d viewM = frustum.ComputeViewMatrix().GetTranspose();
        pxr::GfMatrix4d projM = frustum.ComputeProjectionMatrix().GetTranspose();
        const pxr::UsdPrim cameraPrim = camera.GetPrim();
        pxr::UsdGeomXformCache xformCache(timeCode);
        const pxr::GfMatrix4d cameraLocalToWorld = xformCache.GetLocalToWorldTransform(cameraPrim);
        const pxr::GfVec3d cameraP = cameraLocalToWorld.Transform(pxr::GfVec3d(0.0, 0.0, 0.0));
        const pxr::GfVec3d cameraFwd = cameraLocalToWorld.TransformDir(pxr::GfVec3d(0.0, 0.0, -1.0));
        float3 worldPos = make_float3(float(cameraP[0]), float(cameraP[1]), float(cameraP[2]));
        float3 forward = make_float3(float(cameraFwd[0]), float(cameraFwd[1]), float(cameraFwd[2]));
        const float forwardLen = length(forward);
        if (forwardLen > 1e-8f)
        {
            forward = forward * (1.0f / forwardLen);
        }
        else
        {
            forward = make_float3(0.0f, 0.0f, -1.0f);
        }
        const float verticalFovDegrees =
            float(gfCam.GetFieldOfView(pxr::GfCamera::FOVVertical));
        const pxr::GfRange1d nearFar = frustum.GetNearFar();
        const float nearPlane = float(nearFar.GetMin());

        Camera &uc = scenePool->camera;
        for (int i = 0; i < 4; i++)
        {
            for (int j = 0; j < 4; j++)
            {
                uc.cameraFromWorld.m[i][j] = static_cast<float>(viewM[i][j]);
                uc.clipFromCamera.m[i][j] = static_cast<float>(projM[i][j]);

                printf("%f ", projM[i][j]);
            }
        }

        if (settings)
        {
            pxr::GfVec2i resolution;
            if (settings.GetResolutionAttr().Get(&resolution))
            {
                uc.viewportWidth = resolution[0];
                uc.viewportHeight = resolution[1];
                printf("resolution: %i %i\n", resolution[0], resolution[1]);
            }
        }
        uc.worldPosition = worldPos;
        uc.forward = forward;
        uc.verticalFovDegrees =
            (std::isfinite(verticalFovDegrees) && verticalFovDegrees > 0.0f) ? verticalFovDegrees
                                                                              : 45.0f;
        uc.nearPlane = (std::isfinite(nearPlane) && nearPlane > 0.0f) ? nearPlane : 1.0f;
        uc.hasValidCamera = true;
        uc.path = cameraPrim.GetPath().GetString();

        double shutterOpen, shutterClose;
        camera.GetShutterOpenAttr().Get(&shutterOpen, 0.0);
        camera.GetShutterCloseAttr().Get(&shutterClose, 0.0);
        printf("open: %f close: %f\n", shutterOpen, shutterClose);
    }

    std::unordered_map<std::string, int> materialMap;
    std::vector<pxr::UsdShadeMaterial> materials;

    for (const USDBuildScene &buildScene : buildSceneDAG.scenes)
    {
        for (const USDBuildSceneCurve &curveRef : buildScene.curves)
        {
            const pxr::UsdPrim prim = stage->GetPrimAtPath(pxr::SdfPath(curveRef.path));
            if (!prim || !prim.IsA<pxr::UsdGeomBasisCurves>())
            {
                printf("curve prim missing or invalid: %s\n", curveRef.path.c_str());
                return;
            }
            pxr::UsdShadeMaterial material = GetPrimMaterial(prim);
            AddMaterialToMap(materialMap, materials, material);
        }

        for (const USDBuildSceneMesh &meshRef : buildScene.meshes)
        {
            const pxr::UsdPrim prim = stage->GetPrimAtPath(pxr::SdfPath(meshRef.path));
            if (!prim || !prim.IsA<pxr::UsdGeomMesh>())
            {
                printf("mesh prim missing or invalid: %s\n", meshRef.path.c_str());
                return;
            }

            pxr::UsdGeomMesh mesh(prim);
            const std::vector<pxr::UsdGeomSubset> subsets =
                pxr::UsdGeomSubset::GetAllGeomSubsets(mesh);
            if (subsets.size())
            {
                printf("subsets not handled yet\n");
            }
            else
            {
                pxr::UsdShadeMaterial material = GetPrimMaterial(mesh.GetPrim());
                AddMaterialToMap(materialMap, materials, material);
            }
        }
    }

    printf("num materials: %zi\n", materials.size());

    std::vector<ScenePool::MaterialInfo> materialTextures;
    materialTextures.reserve(materials.size());

    for (pxr::UsdShadeMaterial &material : materials)
    {
        ScenePool::MaterialInfo info = {};
        info.materialPath = material.GetPath().GetString();
        ReadNtcMaterialInfo(material, &info);
        std::vector<std::string> uniqueTexturePaths;

        if (pxr::UsdShadeShader shader = material.ComputeSurfaceSource())
        {
            pxr::TfToken token;
            if (shader.GetShaderId(&token) && token == pxr::TfToken("UsdPreviewSurface"))
            {
                for (const pxr::UsdShadeInput &input : shader.GetInputs())
                {
                    std::string texturePath;
                    std::string swizzle;
                    if (TryGetImageTexturePath(input, texturePath, &swizzle))
                    {
                        AppendUniqueTexturePath(uniqueTexturePaths, texturePath);
                        ScenePool::MaterialTextureInput textureInput = {};
                        textureInput.inputName = input.GetBaseName().GetString();
                        textureInput.texturePath = texturePath;
                        textureInput.swizzle = swizzle;
                        info.textureInputs.push_back(std::move(textureInput));
                    }
                }
            }
            else
            {
                // TODO: support non-UsdPreviewSurface material networks.
                printf("Texture gather: material %s uses unsupported surface shader %s\n",
                       material.GetPath().GetText(),
                       token.GetText());
            }
        }
        else
        {
            printf("Texture gather: material %s has no surface source\n", material.GetPath().GetText());
        }
        materialTextures.push_back(std::move(info));
    }

    for (const ScenePool::MaterialInfo &info : materialTextures)
    {
        printf("material %s image textures: %zu\n",
               info.materialPath.c_str(),
               info.textureInputs.size());
        for (const ScenePool::MaterialTextureInput &textureInput : info.textureInputs)
        {
            printf("  %s (%s)\n", textureInput.texturePath.c_str(), textureInput.inputName.c_str());
        }
        if (!info.ntcDiffuseFile.empty())
        {
            printf("  ntc diffuse: %s (%s)\n",
                   info.ntcDiffuseFile.c_str(),
                   info.ntcDiffuseTextureName.c_str());
        }
    }
    scenePool->materials = std::move(materialTextures);

    int total = 0;
    for (size_t sceneIndex = 0; sceneIndex < buildSceneDAG.scenes.size(); sceneIndex++)
    {
        const USDBuildScene &buildScene = buildSceneDAG.scenes[sceneIndex];
        Scene *outScene = scenePool->scenes[sceneIndex].get();

        for (const USDBuildSceneCurve &curveRef : buildScene.curves)
        {
            const pxr::UsdPrim prim = stage->GetPrimAtPath(pxr::SdfPath(curveRef.path));
            if (!prim || !prim.IsA<pxr::UsdGeomBasisCurves>())
            {
                printf("curve prim missing or invalid: %s\n", curveRef.path.c_str());
                return;
            }

            pxr::UsdGeomBasisCurves curve(prim);
            ProcessUSDBasisCurve(curve, outScene);
            outScene->curves.back().parentFromLocal = curveRef.parentFromLocal;
        }

        outScene->attributes.Reserve(buildScene.meshes.size());
        for (const USDBuildSceneMesh &meshRef : buildScene.meshes)
        {
            const pxr::UsdPrim prim = stage->GetPrimAtPath(pxr::SdfPath(meshRef.path));
            if (!prim || !prim.IsA<pxr::UsdGeomMesh>())
            {
                printf("mesh prim missing or invalid: %s\n", meshRef.path.c_str());
                return;
            }

            pxr::UsdGeomMesh mesh(prim);
            const pxr::UsdShadeMaterial material = GetPrimMaterial(mesh.GetPrim());
            const int materialIndex = GetMaterialIndex(materialMap, material);
            pxr::VtVec3fArray positions;
            pxr::VtIntArray faceIndices;
            pxr::VtIntArray faceCounts;
            pxr::TfToken scheme;

            mesh.GetSubdivisionSchemeAttr().Get(&scheme);
            if (scheme == pxr::UsdGeomTokens->catmullClark)
            {
                ProcessCatmullClarkMesh(mesh, outScene);
                total++;
                continue;
            }
            else if (scheme == pxr::UsdGeomTokens->none)
            {
                ProcessPrimvars(mesh.GetPrim(), 0.0, outScene);
            }
            else
            {
                printf("%s\n", scheme.GetText());
            }

            USD_ASSERT(mesh.GetPointsAttr().Get(&positions, 0.0));
            USD_ASSERT(mesh.GetFaceVertexIndicesAttr().Get(&faceIndices, 0.0));
            USD_ASSERT(mesh.GetFaceVertexCountsAttr().Get(&faceCounts, 0.0));

            int numTriangles = 0;
            for (int faceCount : faceCounts)
            {
                if (faceCount != 3)
                {
                    if (faceCount == 4)
                    {
                        numTriangles += 2;
                    }
                    else
                    {
                        printf("n-gon found\n");
                        YBI_ASSERT(0);
                    }
                }
                else
                {
                    numTriangles++;
                }
            }

            Array<float3> finalPositions(positions);
            Array<int> finalIndices(3 * numTriangles);
            Array<float2> finalTexcoords;
            Array<int> finalTexcoordIndices;

            int inputOffset = 0;
            int finalOffset = 0;
            for (int faceCount : faceCounts)
            {
                if (faceCount == 3)
                {
                    for (int i = 0; i < 3; i++)
                    {
                        int index = faceIndices[inputOffset++];
                        YBI_ASSERT(index < positions.size());
                        finalIndices[finalOffset++] = index;
                    }
                }
                else if (faceCount == 4)
                {
                    int tempIndices[4];
                    for (int i = 0; i < 4; i++)
                    {
                        int index = faceIndices[inputOffset++];
                        YBI_ASSERT(index < positions.size());
                        tempIndices[i] = index;
                    }
                    finalIndices[finalOffset++] = tempIndices[0];
                    finalIndices[finalOffset++] = tempIndices[1];
                    finalIndices[finalOffset++] = tempIndices[2];

                    finalIndices[finalOffset++] = tempIndices[0];
                    finalIndices[finalOffset++] = tempIndices[2];
                    finalIndices[finalOffset++] = tempIndices[3];
                }
                else
                {
                    YBI_ASSERT(0);
                }
            }

            const bool hasTexcoords = BuildTriangulatedMeshTexcoords(
                mesh, faceCounts, faceIndices, numTriangles, &finalTexcoords, &finalTexcoordIndices);

            outScene->meshes.emplace_back(
                std::move(finalPositions), std::move(finalIndices), meshRef.parentFromLocal);
            Mesh &outMesh = outScene->meshes.back();
            outMesh.materialIndex = materialIndex;
            if (hasTexcoords)
            {
                outMesh.texcoords.Resize(finalTexcoords.size());
                if (finalTexcoords.size() > 0)
                {
                    memcpy(outMesh.texcoords.data(),
                           finalTexcoords.data(),
                           sizeof(float2) * finalTexcoords.size());
                }
                outMesh.texcoordIndices.Resize(finalTexcoordIndices.size());
                if (finalTexcoordIndices.size() > 0)
                {
                    memcpy(outMesh.texcoordIndices.data(),
                           finalTexcoordIndices.data(),
                           sizeof(int) * finalTexcoordIndices.size());
                }
            }
        }
    }

    printf("num cat clarks: %i\n", total);
}

YBI_NAMESPACE_END
