#include "io/usd.h"
#include "core/affine_transform.h"
#include "core/allocation.h"
#include "core/arena.h"
#include "core/array.h"
#include "core/assert.h"
#include "core/bit.h"
#include "core/math.h"
#include "core/thread_context.h"
#include "io/usd_helpers.h"
#include "io/usd_material.h"
#include "pxr/base/gf/vec3f.h"
#include "scene/scene.h"
#include <algorithm>
#include <atomic>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <filesystem>
#include <memory>
#include <pxr/base/gf/camera.h>
#include <pxr/base/gf/frustum.h>
#include <pxr/base/tf/stringUtils.h>
#include <pxr/base/vt/array.h>
#include <pxr/base/vt/types.h>
#include <pxr/usd/sdf/assetPath.h>
#include <pxr/usd/usd/attributeQuery.h>
#include <pxr/usd/usd/prim.h>
#include <pxr/usd/usd/primFlags.h>
#include <pxr/usd/usd/primRange.h>
#include <pxr/usd/usd/stage.h>
#include <pxr/usd/usdGeom/basisCurves.h>
#include <pxr/usd/usdGeom/camera.h>
#include <pxr/usd/usdGeom/mesh.h>
#include <pxr/usd/usdGeom/pointInstancer.h>
#include <pxr/usd/usdGeom/primvarsAPI.h>
#include <pxr/usd/usdGeom/subset.h>
#include <pxr/usd/usdGeom/tokens.h>
#include <pxr/usd/usdGeom/xformCache.h>
#include <pxr/usd/usdLux/rectLight.h>
#include <pxr/usd/usdLux/sphereLight.h>
#include <pxr/usd/usdShade/connectableAPI.h>
#include <pxr/usd/usdShade/material.h>
#include <pxr/usd/usdShade/materialBindingAPI.h>
#include <pxr/usd/usdShade/nodeGraph.h>
#include <pxr/usd/usdShade/output.h>
#include <pxr/usd/usdShade/udimUtils.h>
#include <pxr/usd/usdShade/utils.h>
#include <string>
#include <tbb/parallel_for.h>
#include <unordered_set>
#include <utility>
#include <variant>
#include <vector>

namespace ybi
{

static const pxr::TfToken stToken("st", pxr::TfToken::Immortal);
static const pxr::TfToken map1Token("map1", pxr::TfToken::Immortal);
static const pxr::TfToken usdTransform2dToken("UsdTransform2d", pxr::TfToken::Immortal);
static const pxr::TfToken usdPrimvarReaderFloat2Token("UsdPrimvarReader_float2",
                                                      pxr::TfToken::Immortal);
static const pxr::TfToken varnameToken("varname", pxr::TfToken::Immortal);

#define USD_CONCURRENT_MESH_LOAD
#define USD_YBI_ASSERT(expr) expr

constexpr size_t INVALID_INDEX = ~size_t(0);
constexpr uint32_t INVALID_OFFSET = ~0u;

static AttributeType ConvertTokenToAttributeType(const pxr::TfToken &token)
{
    AttributeType type = AttributeType::None;
    if (token == "vertex")
    {
        type = AttributeType::Vertex;
    }
    else if (token == "faceVarying")
    {
        type = AttributeType::FaceVarying;
    }
    else
    {
        YBI_LOGFATAL("Type: %s not supported\n", token.GetText());
    }
    return type;
}

struct USDPrimTraversalData
{
    pxr::UsdPrim prim;
    std::vector<pxr::UsdPrim> meshes;
    std::vector<pxr::UsdPrim> pointInstancers;
    std::vector<pxr::UsdPrim> instances;
    std::vector<pxr::UsdPrim> curves;

    std::vector<pxr::UsdPrim> rectLights;
    std::vector<pxr::UsdPrim> sphereLights;

    std::vector<pxr::UsdGeomCamera> cameras;
};

template <typename T>
MeshBufferInternState<T>::MeshBufferInternState(Arena *arena, size_t numMeshes)
    : map(arena, Max<size_t>(1, NextPowerOfTwo(Max<size_t>(1, numMeshes)))), bufferOffset(0),
      buffers(numMeshes)
{
    indices.Resize(numMeshes);
    for (size_t i = 0; i < indices.size(); i++)
    {
        indices[i] = INVALID_INDEX;
    }
}

template <typename T>
size_t MeshBufferInternState<T>::Intern(Arena *arena, Array<T> &&buffer, size_t bufferIndex)
{
    if (buffer.size() == 0)
    {
        indices[bufferIndex] = INVALID_INDEX;
        return INVALID_INDEX;
    }

    uint64_t hash = MurmurHash64A(
        reinterpret_cast<const uint8_t *>(buffer.data()), buffer.size() * sizeof(T), 0);
    size_t duplicateIndex = map.GetOrAdd(
        arena,
        bufferIndex,
        hash,
        [&](size_t index) {
            YBI_ERROR(index < buffers.size(), "%zu %zu\n", index, buffers.size());
            const InternedMeshBuffer<T> &other = buffers[index];

            if (hash != other.hash || other.buffer.size() != buffer.size())
                return false;

            return memcmp(buffer.data(), other.buffer.data(), sizeof(T) * buffer.size()) == 0;
        },
        [&](size_t index) {
            InternedMeshBuffer<T> &entry = buffers[index];
            entry.hash = hash;
            entry.buffer = std::move(buffer);
            bufferOffset.fetch_add(1, std::memory_order_relaxed);
        });

    indices[bufferIndex] = duplicateIndex;
    return duplicateIndex;
}

template <typename T>
void MeshBufferInternState<T>::MoveSortedBuffersTo(std::vector<Array<T>> *dst)
{
    // NOTE: during concurrent loads, sort the buffers to maintain a deterministic order
    // across separate program runs
    struct SortEntry
    {
        uint64_t hash;
        size_t oldIndex;
        Array<T> buffer;
    };

    std::vector<SortEntry> sortedBuffers;
    sortedBuffers.reserve(bufferOffset.load(std::memory_order_relaxed));

    oldToNew.Resize(buffers.size());
    for (int i = 0; i < buffers.size(); i++)
    {
        oldToNew[i] = INVALID_OFFSET;
    }

    for (size_t oldIndex = 0; oldIndex < buffers.size(); oldIndex++)
    {
        InternedMeshBuffer<T> &entry = buffers[oldIndex];
        if (entry.buffer.size())
        {
            sortedBuffers.push_back({entry.hash, oldIndex, std::move(entry.buffer)});
        }
    }

    std::sort(
        sortedBuffers.begin(), sortedBuffers.end(), [](const SortEntry &a, const SortEntry &b) {
            if (a.hash != b.hash)
                return a.hash < b.hash;
            return a.oldIndex < b.oldIndex;
        });

    dst->clear();
    dst->reserve(sortedBuffers.size());
    for (uint32_t newIndex = 0; newIndex < sortedBuffers.size(); newIndex++)
    {
        SortEntry &entry = sortedBuffers[newIndex];
        oldToNew[entry.oldIndex] = newIndex;
        dst->push_back(std::move(entry.buffer));
    }
}

template <typename T>
uint32_t MeshBufferInternState<T>::ResolveIndex(size_t index) const
{
    if (index == INVALID_INDEX || index == INVALID_OFFSET)
    {
        return INVALID_OFFSET;
    }

    uint32_t newIndex = oldToNew[index];
    YBI_ASSERT(newIndex != INVALID_OFFSET);
    return newIndex;
}

MeshBufferInternStates::MeshBufferInternStates(Arena *arena, size_t numMeshes)
    : positions(arena, numMeshes), normals(arena, numMeshes), uvs(arena, numMeshes),
      uvIndices(arena, numMeshes), indices(arena, numMeshes), faceMaterialIndices(arena, numMeshes)
{
}

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

static std::string DebugMaterialValue(const pxr::UsdShadeInput &input)
{
    pxr::VtValue value;
    if (!input.Get(&value))
    {
        return {};
    }

    std::string result = pxr::TfStringify(value);
    constexpr size_t maxLength = 160;
    if (result.size() > maxLength)
    {
        result.resize(maxLength);
        result += "...";
    }
    return result;
}

static void DumpMaterialShader(const pxr::UsdShadeShader &shader,
                               int depth,
                               std::unordered_set<std::string> &visitedShaders)
{
    if (!shader)
    {
        return;
    }

    const std::string shaderPath = shader.GetPath().GetString();
    pxr::TfToken shaderId;
    shader.GetShaderId(&shaderId);

    printf("%*sShader %s [%s]", depth * 2, "", shaderPath.c_str(), shaderId.GetText());
    if (!visitedShaders.insert(shaderPath).second)
    {
        printf(" (already shown)\n");
        return;
    }
    printf("\n");

    for (const pxr::UsdShadeInput &input : shader.GetInputs())
    {
        const std::string value = DebugMaterialValue(input);
        const pxr::UsdShadeSourceInfoVector sources = input.GetConnectedSources();

        printf("%*s+- input %s : %s",
               (depth + 1) * 2,
               "",
               input.GetBaseName().GetText(),
               input.GetTypeName().GetAsToken().GetText());
        if (!value.empty())
        {
            printf(" = %s", value.c_str());
        }

        if (sources.empty())
        {
            printf("\n");
            continue;
        }
        printf("\n");

        for (const pxr::UsdShadeConnectionSourceInfo &source : sources)
        {
            const pxr::UsdPrim sourcePrim = source.source.GetPrim();
            printf("%*s<- %s.%s (%s)\n",
                   (depth + 2) * 2,
                   "",
                   sourcePrim.GetPath().GetText(),
                   source.sourceName.GetText(),
                   source.sourceType == pxr::UsdShadeAttributeType::Output ? "output" : "input");

            if (sourcePrim.IsA<pxr::UsdShadeShader>())
            {
                DumpMaterialShader(pxr::UsdShadeShader(sourcePrim), depth + 3, visitedShaders);
                continue;
            }

            if (sourcePrim.IsA<pxr::UsdShadeNodeGraph>())
            {
                const pxr::UsdShadeNodeGraph nodeGraph(sourcePrim);
                const pxr::UsdShadeOutput output = nodeGraph.GetOutput(source.sourceName);
                for (const pxr::UsdAttribute &attribute :
                     pxr::UsdShadeUtils::GetValueProducingAttributes(output))
                {
                    const pxr::UsdPrim producerPrim = attribute.GetPrim();
                    if (producerPrim.IsA<pxr::UsdShadeShader>())
                    {
                        DumpMaterialShader(
                            pxr::UsdShadeShader(producerPrim), depth + 3, visitedShaders);
                    }
                }
                continue;
            }

            if (source.sourceType == pxr::UsdShadeAttributeType::Input)
            {
                const pxr::UsdShadeInput sourceInput = source.source.GetInput(source.sourceName);
                const std::string sourceValue = DebugMaterialValue(sourceInput);
                if (!sourceValue.empty())
                {
                    printf("%*svalue = %s\n", (depth + 3) * 2, "", sourceValue.c_str());
                }
            }
        }
    }
}

static void DumpMaterialGraph(const pxr::UsdShadeMaterial &material,
                              const pxr::UsdShadeShader &surfaceShader)
{
    printf("\n========== USD MATERIAL GRAPH ==========\n");
    printf("Material %s\n", material.GetPath().GetText());

    std::unordered_set<std::string> visitedShaders;
    DumpMaterialShader(surfaceShader, 1, visitedShaders);

    printf("========================================\n\n");
}

static uint32_t AddMaterialToMap(std::unordered_map<std::string, size_t> &materialMap,
                                 std::vector<pxr::UsdShadeMaterial> &materials,
                                 const pxr::UsdShadeMaterial &material)
{
    std::string materialPath = material.GetPath().GetString();
    auto it = materialMap.find(materialPath);
    if (it == materialMap.end())
    {
        materials.push_back(material);
        materialMap.insert({materialPath, materials.size() - 1});
        return static_cast<uint32_t>(materials.size() - 1);
    }
    else
    {
        return static_cast<uint32_t>(it->second);
    }
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

template <typename OutputT, typename Values, typename Convert>
static Array<OutputT> TriangulateFaceVaryingValues(const Values &values,
                                                   const pxr::VtIntArray &faceCounts,
                                                   Convert &&convert)
{
    Array<OutputT> triangles;
    size_t faceVertexOffset = 0;

    for (const int faceCount : faceCounts)
    {
        YBI_ASSERT(faceCount >= 3);
        YBI_ASSERT(faceVertexOffset + static_cast<size_t>(faceCount) <= values.size());

        for (size_t i = 1; i < static_cast<size_t>(faceCount) - 1; ++i)
        {
            triangles.EmplaceBack(convert(values[faceVertexOffset]));
            triangles.EmplaceBack(convert(values[faceVertexOffset + i]));
            triangles.EmplaceBack(convert(values[faceVertexOffset + i + 1]));
        }

        faceVertexOffset += static_cast<size_t>(faceCount);
    }

    YBI_ASSERT(faceVertexOffset == values.size());
    return triangles;
}

template <typename Values>
static Array<typename Values::value_type>
TriangulateFaceVaryingValues(const Values &values, const pxr::VtIntArray &faceCounts)
{
    return TriangulateFaceVaryingValues<typename Values::value_type>(
        values, faceCounts, [](const auto &value) { return value; });
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

AttributeType USDSceneLoader::ProcessPrimvars(const pxr::UsdGeomMesh &mesh,
                                              pxr::VtArray<pxr::GfVec2f> *outUVs,
                                              pxr::VtIntArray *outUVIndices,
                                              pxr::UsdTimeCode timeCode) const
{
    outUVs->clear();
    outUVIndices->clear();

    pxr::UsdGeomPrimvarsAPI primvars(mesh);

    pxr::UsdGeomPrimvar stPrimvar = primvars.GetPrimvar(stToken);

    if (!stPrimvar)
    {
        return AttributeType::None;
    }

    pxr::VtVec2fArray stValues;
    if (!stPrimvar.Get(&stValues, timeCode))
    {
        return AttributeType::None;
    }
    if (stValues.empty())
    {
        return AttributeType::None;
    }

    pxr::VtIntArray indices;
    const bool hasIndices =
        stPrimvar.IsIndexed() && stPrimvar.GetIndices(&indices, timeCode) && !indices.empty();

    pxr::TfToken interpolationToken = stPrimvar.GetInterpolation();

    AttributeType type = ConvertTokenToAttributeType(interpolationToken);

    *outUVs = stValues;
    if (hasIndices)
    {
        *outUVIndices = indices;
    }

    return type;
}

static bool HasAllZeroNormals(const pxr::VtVec3fArray &normals)
{
    for (const pxr::GfVec3f &normal : normals)
    {
        const float lengthSquared =
            normal[0] * normal[0] + normal[1] * normal[1] + normal[2] * normal[2];
        if (lengthSquared > 1e-20f)
        {
            return false;
        }
    }

    return true;
}

AttributeType USDSceneLoader::ProcessNormals(const pxr::UsdGeomMesh &mesh,
                                             pxr::VtVec3fArray *outNormals,
                                             pxr::UsdTimeCode timeCode) const
{
    pxr::TfToken interpolationToken;

    pxr::UsdGeomPrimvarsAPI primvars(mesh);
    pxr::UsdGeomPrimvar normalsPrimvar = primvars.GetPrimvar(pxr::UsdGeomTokens->normals);

    if (normalsPrimvar)
    {
        pxr::VtValue values;
        if (!normalsPrimvar.ComputeFlattened(&values, timeCode))
        {
            return AttributeType::None;
        }
        if (!values.IsHolding<pxr::VtVec3fArray>())
        {
            YBI_LOGFATAL("primvars:normals are not GfVec3f\n");
            return AttributeType::None;
        }

        *outNormals = values.Get<pxr::VtVec3fArray>();
        if (outNormals->empty() || HasAllZeroNormals(*outNormals))
        {
            return AttributeType::None;
        }

        interpolationToken = normalsPrimvar.GetInterpolation();
    }

    if (interpolationToken.IsEmpty() && mesh.GetNormalsAttr().Get(outNormals, timeCode))
    {
        if (!outNormals->empty() && !HasAllZeroNormals(*outNormals))
        {
            interpolationToken = mesh.GetNormalsInterpolation();
        }
    }

    if (interpolationToken.IsEmpty())
    {
        return AttributeType::None;
    }

    if (interpolationToken == pxr::UsdGeomTokens->vertex)
    {
        return AttributeType::Vertex;
    }
    else if (interpolationToken == pxr::UsdGeomTokens->faceVarying)
    {
        return AttributeType::FaceVarying;
    }
    else
    {
        YBI_LOGFATAL("Need to handle other normal interpolation types\n");
        return AttributeType::None;
    }
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

    pxr::VtVec2fArray usdUvs;
    pxr::VtIntArray usdUvIndices;
    Array<Vec2f> uvs;
    Array<int> uvIndices;
    AttributeType uvType = ProcessPrimvars(usdMesh, &usdUvs, &usdUvIndices, timeCode);
    uvs.Reserve(usdUvs.size());
    switch (uvType)
    {
        case AttributeType::FaceVarying:
        {
            if (usdUvIndices.size())
            {
                YBI_ASSERT(usdUvIndices.size() == faceIndices.size());
                uvs.Resize(usdUvs.size());
                memcpy(uvs.data(), usdUvs.data(), sizeof(Vec2f) * usdUvs.size());
                uvIndices = TriangulateFaceVaryingValues(usdUvIndices, faceCounts);
            }
            else
            {
                YBI_ASSERT(usdUvs.size() == faceIndices.size());
                uvs = TriangulateFaceVaryingValues<Vec2f>(
                    usdUvs, faceCounts, [](const pxr::GfVec2f &value) {
                        return Vec2f(value[0], value[1]);
                    });
            }
            break;
        }
        case AttributeType::Vertex:
        {
            uvs.Resize(usdUvs.size());
            memcpy(uvs.data(), usdUvs.data(), sizeof(Vec2f) * usdUvs.size());
            if (usdUvIndices.size())
            {
                YBI_ASSERT(usdUvIndices.size() == usdPositions.size());
                uvIndices.Resize(usdUvIndices.size());
                memcpy(uvIndices.data(), usdUvIndices.data(), sizeof(int) * usdUvIndices.size());
            }
            else
            {
                YBI_ASSERT(usdUvs.size() == usdPositions.size());
            }
            break;
        }
        case AttributeType::None:
        {
            break;
        }
        default:
        {
            YBI_LOGFATAL("Only vertex/face varying UVs supported currently.");
            break;
        }
    }

    pxr::VtVec3fArray usdNormals;
    Array<Vec3f> normals;
    AttributeType normalType = ProcessNormals(usdMesh, &usdNormals, timeCode);
    normals.Reserve(usdNormals.size());
    switch (normalType)
    {
        case AttributeType::FaceVarying:
        {
            YBI_ASSERT(usdNormals.size() == faceIndices.size());
            normals = TriangulateFaceVaryingValues<Vec3f>(
                usdNormals, faceCounts, [](const pxr::GfVec3f &value) {
                    return Vec3f(value[0], value[1], value[2]);
                });
            break;
        }
        case AttributeType::Vertex:
        {
            YBI_ASSERT(usdNormals.size() == usdPositions.size());
            normals.Resize(usdNormals.size());
            memcpy(normals.data(), usdNormals.data(), sizeof(Vec3f) * usdNormals.size());

            break;
        }
        case AttributeType::None:
        {
            break;
        }
        default:
        {
            YBI_LOGFATAL("Only vertex/face varying normals supported currently.");
            break;
        }
    }

    Array<Vec3f> positions(usdPositions.size());
    memcpy(positions.data(), usdPositions.data(), sizeof(Vec3f) * usdPositions.size());
    Array<int> triangleIndices = TriangulateFaceVaryingValues(faceIndices, faceCounts);

    printf("%zu %zu %zu\n", positions.size(), normals.size(), uvs.size());
    size_t index = meshOffset + meshIndex;
    BufferHandle positionHandle{static_cast<uint32_t>(scene->container->meshPositions.size())};
    scene->container->meshPositions.push_back(std::move(positions));

    BufferHandle normalHandle;
    if (normals.size())
    {
        normalHandle.index = static_cast<uint32_t>(scene->container->meshNormals.size());
        scene->container->meshNormals.push_back(std::move(normals));
    }

    BufferHandle uvHandle;
    if (uvs.size())
    {
        uvHandle.index = static_cast<uint32_t>(scene->container->meshUVs.size());
        scene->container->meshUVs.push_back(std::move(uvs));
    }

    BufferHandle uvIndexHandle;
    if (uvIndices.size())
    {
        uvIndexHandle.index = static_cast<uint32_t>(scene->container->meshUVIndices.size());
        scene->container->meshUVIndices.push_back(std::move(uvIndices));
    }

    BufferHandle indexHandle{static_cast<uint32_t>(scene->container->meshTriangleIndices.size())};
    scene->container->meshTriangleIndices.push_back(std::move(triangleIndices));

    Mesh mesh(
        positionHandle, normalHandle, normalType, uvHandle, uvIndexHandle, uvType, indexHandle);
    scene->container->meshes[index] = std::move(mesh);

    bool resetsXformStack = false;
    pxr::GfMatrix4d localFromMesh =
        cache.ComputeRelativeTransform(prim, sceneRoot, &resetsXformStack);
    Instance &instance = scene->instances[meshIndex];
    instance.transform = ConvertAffineTransform(localFromMesh);
    instance.handle.type = GeometryType::Mesh;
    instance.handle.index = index;
}

// NOTE: Meshes may have duplicated buffers (positions, indices, etc). This minimizes
// final mesh allocation size.
void USDSceneLoader::ProcessPolygonMeshConcurrent(Arena *arena,
                                                  Scene *scene,
                                                  pxr::UsdPrim sceneRoot,
                                                  size_t meshOffset,
                                                  size_t meshIndex,
                                                  const pxr::UsdPrim prim,
                                                  MeshBufferInternStates *internStates,
                                                  pxr::UsdGeomXformCache &cache,
                                                  pxr::UsdTimeCode timeCode)
{
    pxr::UsdGeomMesh usdMesh(prim);

    pxr::UsdGeomImageable imageable(prim);

    pxr::VtVec3fArray usdPositions;
    pxr::VtIntArray faceIndices;
    pxr::VtIntArray faceCounts;
    pxr::TfToken subdivisionScheme;

    pxr::UsdAttribute pointsAttr = usdMesh.GetPointsAttr();
    pxr::UsdAttributeQuery pointsQuery(pointsAttr);
    bool isMaybeTimeVarying = pointsQuery.ValueMightBeTimeVarying();

    bool result = pointsAttr.Get(&usdPositions, timeCode);
    YBI_ASSERT(result);
    YBI_ASSERT(usdPositions.size());

#if 0
    if (isMaybeTimeVarying)
    {
        std::vector<double> times;
        bool noError = pointsAttr.GetTimeSamples(&times);
        printf("num samples: %zu\n", times.size());
        if (noError)
        {
            for (double time : times)
            {
                int maxDeltaXorA = 0;
                int maxDeltaXorB = 0;
                int maxDeltaXorC = 0;

                uint32_t maxDeltaDeltaA = 0;
                uint32_t maxDeltaDeltaB = 0;
                uint32_t maxDeltaDeltaC = 0;

                printf("time: %f\n", time);
                pxr::VtVec3fArray usdTimePositions;
                result = pointsAttr.Get(&usdTimePositions, pxr::UsdTimeCode(time));
                YBI_ASSERT(result);

                YBI_ASSERT(usdTimePositions.size() == usdPositions.size());
                for (int i = 0; i < usdTimePositions.size(); i++)
                {
                    pxr::GfVec3f p0 = usdPositions[i];
                    pxr::GfVec3f p1 = usdTimePositions[i];

                    pxr::GfVec3f test = p1 - p0;

                    int a = FloatToBits(p0[0]);
                    int b = FloatToBits(p0[1]);
                    int c = FloatToBits(p0[2]);

                    int d = FloatToBits(p1[0]);
                    int e = FloatToBits(p1[1]);
                    int f = FloatToBits(p1[2]);

                    // printf("%f %f %f\n", test[0], test[1], test[2]);
                    // printf("%u %u %u\n", a ^ d, b ^ e, c ^ f);
                    // printf("%u %u %u\n",
                    //        ZigZagEncode(d - a),
                    //        ZigZagEncode(e - b),
                    //        ZigZagEncode(f - c));

                    maxDeltaXorA = Max(maxDeltaXorA, a ^ d);
                    maxDeltaXorB = Max(maxDeltaXorB, b ^ e);
                    maxDeltaXorC = Max(maxDeltaXorC, c ^ f);

                    maxDeltaDeltaA = Max(maxDeltaDeltaA, ZigZagEncode(d - a));
                    maxDeltaDeltaB = Max(maxDeltaDeltaB, ZigZagEncode(e - b));
                    maxDeltaDeltaC = Max(maxDeltaDeltaC, ZigZagEncode(f - c));
                }
                printf("xor: %i %i %i, delta: %u %u %u\n",
                       maxDeltaXorA,
                       maxDeltaXorB,
                       maxDeltaXorC,
                       maxDeltaDeltaA,
                       maxDeltaDeltaB,
                       maxDeltaDeltaC);
            }
        }
        printf("end\n");
    }
#endif

    Array<Vec3f> positions(usdPositions.size());
    memcpy(positions.data(), usdPositions.data(), sizeof(Vec3f) * usdPositions.size());

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

    pxr::VtVec2fArray usdUvs;
    pxr::VtIntArray usdUvIndices;
    Array<Vec2f> uvs;
    Array<int> uvIndices;
    AttributeType uvType = ProcessPrimvars(usdMesh, &usdUvs, &usdUvIndices, timeCode);
    uvs.Reserve(usdUvs.size());
    switch (uvType)
    {
        case AttributeType::FaceVarying:
        {
            if (usdUvIndices.size())
            {
                YBI_ASSERT(usdUvIndices.size() == faceIndices.size());
                uvs.Resize(usdUvs.size());
                memcpy(uvs.data(), usdUvs.data(), sizeof(Vec2f) * usdUvs.size());
                uvIndices = TriangulateFaceVaryingValues(usdUvIndices, faceCounts);
            }
            else
            {
                YBI_ASSERT(usdUvs.size() == faceIndices.size());
                uvs = TriangulateFaceVaryingValues<Vec2f>(
                    usdUvs, faceCounts, [](const pxr::GfVec2f &value) {
                        return Vec2f(value[0], value[1]);
                    });
            }
            break;
        }
        case AttributeType::Vertex:
        {
            uvs.Resize(usdUvs.size());
            memcpy(uvs.data(), usdUvs.data(), sizeof(Vec2f) * usdUvs.size());
            if (usdUvIndices.size())
            {
                YBI_ASSERT(usdUvIndices.size() == usdPositions.size());
                uvIndices.Resize(usdUvIndices.size());
                memcpy(uvIndices.data(), usdUvIndices.data(), sizeof(int) * usdUvIndices.size());
            }
            else
            {
                YBI_ASSERT(usdUvs.size() == usdPositions.size());
            }
            break;
        }
        case AttributeType::None:
        {
            break;
        }
        default:
        {
            YBI_LOGFATAL("Only vertex/face varying UVs supported currently.");
            break;
        }
    }

    pxr::VtVec3fArray usdNormals;
    Array<Vec3f> normals;
    AttributeType normalType = ProcessNormals(usdMesh, &usdNormals, timeCode);
    normals.Reserve(usdNormals.size());
    switch (normalType)
    {
        case AttributeType::FaceVarying:
        {
            YBI_ASSERT(usdNormals.size() == faceIndices.size());
            normals = TriangulateFaceVaryingValues<Vec3f>(
                usdNormals, faceCounts, [](const pxr::GfVec3f &value) {
                    return Vec3f(value[0], value[1], value[2]);
                });
            break;
        }
        case AttributeType::Vertex:
        {
            YBI_ASSERT(usdNormals.size() == usdPositions.size());
            normals.Resize(usdNormals.size());
            memcpy(normals.data(), usdNormals.data(), sizeof(Vec3f) * usdNormals.size());

            break;
        }
        case AttributeType::None:
        {
            break;
        }
        default:
        {
            YBI_LOGFATAL("Only vertex/face varying normals supported currently.");
            break;
        }
    }

    Array<int> triangleIndices = TriangulateFaceVaryingValues(faceIndices, faceCounts);

    size_t index = meshOffset + meshIndex;

    BufferHandle positionHandle{
        static_cast<uint32_t>(internStates->positions.Intern(arena, std::move(positions), index))};
    BufferHandle normalHandle{
        static_cast<uint32_t>(internStates->normals.Intern(arena, std::move(normals), index))};
    BufferHandle uvHandle{
        static_cast<uint32_t>(internStates->uvs.Intern(arena, std::move(uvs), index))};
    BufferHandle uvIndexHandle{
        static_cast<uint32_t>(internStates->uvIndices.Intern(arena, std::move(uvIndices), index))};
    BufferHandle indexHandle{static_cast<uint32_t>(
        internStates->indices.Intern(arena, std::move(triangleIndices), index))};

    Mesh mesh(
        positionHandle, normalHandle, normalType, uvHandle, uvIndexHandle, uvType, indexHandle);
    scene->container->meshes[index] = std::move(mesh);

    bool resetsXformStack = false;
    pxr::GfMatrix4d localFromMesh =
        cache.ComputeRelativeTransform(prim, sceneRoot, &resetsXformStack);
    Instance &instance = scene->instances[meshIndex];
    instance.transform = ConvertAffineTransform(localFromMesh);
    instance.handle.type = GeometryType::Mesh;
    instance.handle.index = index;
}

void USDSceneLoader::GetFaceMaterialIndices(Array<int> &faceMaterialIndices,
                                            const pxr::VtIntArray &faceCounts,
                                            const std::vector<pxr::UsdGeomSubset> &subsets,
                                            std::unordered_map<std::string, size_t> &materialMap,
                                            std::vector<pxr::UsdShadeMaterial> &materials,
                                            pxr::UsdTimeCode timeCode) const
{

    ScratchArena faceOffsetScratch;
    View<int> faceOffsets = faceOffsetScratch.arena->PushArray<int>(faceCounts.size() + 1);
    int faceOffset = 0;

    for (int faceIndex = 0; faceIndex < faceCounts.size(); faceIndex++)
    {
        faceOffsets[faceIndex] = faceOffset;
        faceOffset += faceCounts[faceIndex] - 2;
    }
    faceOffsets[faceCounts.size()] = faceOffset;
    faceMaterialIndices.Resize(faceOffset);
    for (int &index : faceMaterialIndices)
    {
        index = -1;
    }

    for (pxr::UsdGeomSubset subset : subsets)
    {
        pxr::VtIntArray subsetFaces;
        subset.GetIndicesAttr().Get(&subsetFaces, timeCode);

        pxr::TfToken elementType;
        subset.GetElementTypeAttr().Get(&elementType, timeCode);
        YBI_ERROR(elementType == pxr::UsdGeomTokens->face,
                  "Subset element type %s not supported\n",
                  elementType.GetString().c_str());

        pxr::UsdShadeMaterial material = GetPrimMaterial(subset.GetPrim());

        if (!material)
        {
            // printf("%s is missing a material\n", prim.GetPath().GetString().c_str());
            continue;
        }

        uint32_t materialIndex = AddMaterialToMap(materialMap, materials, material);

        for (int face : subsetFaces)
        {
            for (int offset = faceOffsets[face]; offset < faceOffsets[face + 1]; offset++)
            {
                YBI_ASSERT(faceMaterialIndices[offset] == -1 ||
                           faceMaterialIndices[offset] == materialIndex);
                faceMaterialIndices[offset] = static_cast<int>(materialIndex);
            }
        }
    }
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
        else if (prim.IsA<pxr::UsdLuxRectLight>())
        {
            data->rectLights.push_back(prim);
        }
        else if (prim.IsA<pxr::UsdLuxSphereLight>())
        {
            data->sphereLights.push_back(prim);
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

static void ParseSTInput(const pxr::UsdShadeShader &shader)
{
    pxr::UsdShadeInput stInput = shader.GetInput(stToken);
    if (!stInput)
    {
        for (pxr::UsdShadeInput input : shader.GetInputs())
        {
            printf("%s\n", input.GetFullName().GetText());
        }
        return;
    }

    pxr::UsdShadeConnectableAPI sourceNode;
    pxr::TfToken sourceOutputName;
    pxr::UsdShadeAttributeType sourceType;

    if (stInput.GetConnectedSource(&sourceNode, &sourceOutputName, &sourceType))
    {
        pxr::UsdPrim sourcePrim = sourceNode.GetPrim();
        pxr::UsdShadeShader connectedShader(sourcePrim);

        pxr::TfToken idToken;
        connectedShader.GetShaderId(&idToken);

        if (idToken == usdPrimvarReaderFloat2Token)
        {
            pxr::UsdShadeInput varnameInput = connectedShader.GetInput(varnameToken);

            pxr::TfToken varname;
            varnameInput.Get(&varname);

            // for (auto &source : connectedShader.GetInput(varnameToken).GetConnectedSources())
            // {
            //     printf("boom: %s\n", source.sourceName.GetText());
            //     // pxr::UsdShadeShader s = pxr::UsdShadeShader(source.source.GetPrim());
            //     // printf("boom %s\n", s.GetInput(source.sourceName);
            // }
        }
        if (idToken == usdTransform2dToken)
        {
            printf("woohoo\n");
        }
    }
}

static PrimvarValue LoadPrimvarValue(const pxr::VtValue &value)
{
    if (value.IsHolding<bool>())
    {
        return value.UncheckedGet<bool>();
    }
    if (value.IsHolding<int>())
    {
        return value.UncheckedGet<int>();
    }
    if (value.IsHolding<float>())
    {
        return value.UncheckedGet<float>();
    }
    if (value.IsHolding<pxr::VtVec2fArray>())
    {
        const pxr::VtVec2fArray &source = value.UncheckedGet<pxr::VtVec2fArray>();
        Array<Vec2f> result(source.size());
        for (size_t i = 0; i < source.size(); ++i)
            result[i] = Vec2f(source[i][0], source[i][1]);
        return result;
    }
    if (value.IsHolding<pxr::VtVec3fArray>())
    {
        const pxr::VtVec3fArray &source = value.UncheckedGet<pxr::VtVec3fArray>();
        Array<Vec3f> result(source.size());
        for (size_t i = 0; i < source.size(); ++i)
            result[i] = Vec3f(source[i][0], source[i][1], source[i][2]);
        return result;
    }
    if (value.IsHolding<pxr::VtVec4fArray>())
    {
        const pxr::VtVec4fArray &source = value.UncheckedGet<pxr::VtVec4fArray>();
        Array<Vec4f> result(source.size());
        for (size_t i = 0; i < source.size(); ++i)
            result[i] = Vec4f(source[i][0], source[i][1], source[i][2], source[i][3]);
        return result;
    }

    YBI_LOGFATAL("Unsupported primvar value type: %s\n", value.GetTypeName().c_str());
    return {};
}

static bool LoadPrimvar(const pxr::UsdStageRefPtr &stage,
                        const pxr::SdfPath &primPath,
                        pxr::UsdTimeCode timeCode,
                        const std::string &name,
                        MeshAttribute &outAttribute)
{
    const pxr::UsdPrim prim = stage->GetPrimAtPath(primPath);
    if (!prim)
    {
        YBI_LOGFATAL("Cannot load geometry attribute '%s': source prim '%s' no longer exists.\n",
                     name.c_str(),
                     primPath.GetText());
        return false;
    }

    pxr::UsdGeomPrimvar primvar =
        pxr::UsdGeomPrimvarsAPI(prim).FindPrimvarWithInheritance(pxr::TfToken(name));

    if (!primvar)
    {
        YBI_LOGFATAL("Material graph references primvar '%s', but geometry prim '%s' does not "
                     "define it.\n",
                     name.c_str(),
                     prim.GetPath().GetText());
        return false;
    }

    pxr::VtValue value;
    if (!primvar.Get(&value, timeCode))
    {
        YBI_LOGFATAL("Material graph references primvar '%s', but geometry prim '%s' has no "
                     "value for it.\n",
                     name.c_str(),
                     prim.GetPath().GetText());
        return false;
    }

    pxr::TfToken interpolationToken = primvar.GetInterpolation();
    outAttribute.type = ConvertTokenToAttributeType(interpolationToken);
    outAttribute.value = LoadPrimvarValue(value);

    pxr::VtIntArray indices;
    const bool hasIndices =
        primvar.IsIndexed() && primvar.GetIndices(&indices, timeCode) && !indices.empty();

    if (outAttribute.type == AttributeType::FaceVarying)
    {
        pxr::UsdGeomMesh mesh(prim);
        pxr::VtIntArray faceCounts;
        YBI_ERROR(mesh.GetFaceVertexCountsAttr().Get(&faceCounts, timeCode),
                  "Could not read face vertex counts for primvar '%s' on mesh '%s'.\n",
                  name.c_str(),
                  prim.GetPath().GetText());

        if (hasIndices)
        {
            // Indexed primvars keep their value table intact; their per-corner indices must
            // follow the triangulated corner order.
            Array<int> sourceIndices(indices.size());
            memcpy(sourceIndices.data(), indices.data(), indices.size() * sizeof(int));
            outAttribute.indices = TriangulateFaceVaryingValues(sourceIndices, faceCounts);
        }
        else
        {
            std::visit(
                [&faceCounts](auto &values) {
                    using T = std::decay_t<decltype(values)>;
                    if constexpr (std::is_same_v<T, Array<Vec2f>> ||
                                  std::is_same_v<T, Array<Vec3f>> ||
                                  std::is_same_v<T, Array<Vec4f>>)
                    {
                        values = TriangulateFaceVaryingValues(values, faceCounts);
                    }
                    else
                    {
                        YBI_LOGFATAL("Face-varying primvars must have an array value.\n");
                    }
                },
                outAttribute.value);
        }
    }
    else if (hasIndices)
    {
        outAttribute.indices.Resize(indices.size());
        memcpy(outAttribute.indices.data(),
               indices.data(),
               outAttribute.indices.size() * sizeof(int));
    }

    return outAttribute.type != AttributeType::None;
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

Scene *USDSceneLoader::LoadScene(Arena *arena,
                                 DataContainer *container,
                                 SceneImportContext *importContext)
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

    printf("camera name: %s\n", camera.GetPath().GetText());
    pxr::UsdPrim cameraPrim = camera.GetPrim();
    pxr::GfCamera gfCam = camera.GetCamera(timeCode);
    const float verticalFovRadians =
        DegreesToRadians(float(gfCam.GetFieldOfView(pxr::GfCamera::FOVVertical)));
    pxr::GfFrustum frustum = gfCam.GetFrustum();
    pxr::GfMatrix4d viewMatrix = frustum.ComputeViewMatrix();
    pxr::GfMatrix4d projMatrix = frustum.ComputeProjectionMatrix();

    printf("fov: %f %f\n",
           verticalFovRadians,
           float(gfCam.GetFieldOfView(pxr::GfCamera::FOVVertical)));

    container->camera.cameraFromWorld = ConvertAffineTransform(viewMatrix);
    AffineTransform renderFromCamera;
    bool inverseResult = container->camera.cameraFromWorld.TryInverse(&renderFromCamera);
    YBI_ASSERT(inverseResult);
    container->camera.position = renderFromCamera.GetTranslation();
    container->camera.verticalFovRadians = verticalFovRadians;
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

    std::unordered_map<std::string, size_t> materialMap;
    std::vector<size_t> sceneMeshOffsets(traversalDatas.size());

    // Set up concurrent hash map for deduplicated buffers while loading
    ScratchArena scratch;
    size_t totalNumMeshes = 0;
    for (size_t i = 0; i < traversalDatas.size(); i++)
    {
        totalNumMeshes += traversalDatas[i].meshes.size();
    }

#ifdef USD_CONCURRENT_MESH_LOAD
    tbb::enumerable_thread_specific<Arena> threadArenas([&]() { return Arena(); });
    MeshBufferInternStates internStates(scratch.arena, totalNumMeshes);
#endif

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
        sceneMeshOffsets[i] = meshOffset;
        size_t dstInstanceIndex = data.meshes.size();
        scene->instances.Resize(totalInstanceCount + data.meshes.size());
        scene->container->meshes.resize(meshOffset + data.meshes.size());

        tbb::parallel_for(tbb::blocked_range<size_t>(0, data.meshes.size()),
                          [&](const tbb::blocked_range<size_t> &range) {
                              pxr::UsdGeomXformCache &cache = xformCaches.local();
                              for (size_t i = range.begin(); i != range.end(); i++)
                              {
                                  pxr::UsdPrim prim = data.meshes[i];

#ifdef USD_CONCURRENT_MESH_LOAD
                                  ProcessPolygonMeshConcurrent(&threadArenas.local(),
                                                               scene,
                                                               data.prim,
                                                               meshOffset,
                                                               i,
                                                               prim,
                                                               &internStates,
                                                               cache,
                                                               timeCode);
#else
                                  // TODO IMPORTANT: this doesn't work multithreaded anymore
                                  ProcessPolygonMesh(scene, data.prim, meshOffset, i, prim, cache, timeCode);
#endif
                              }
                          });

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

#ifdef USD_CONCURRENT_MESH_LOAD
    internStates.positions.MoveSortedBuffersTo(&container->meshPositions);
    internStates.normals.MoveSortedBuffersTo(&container->meshNormals);
    internStates.uvs.MoveSortedBuffersTo(&container->meshUVs);
    internStates.uvIndices.MoveSortedBuffersTo(&container->meshUVIndices);
    internStates.indices.MoveSortedBuffersTo(&container->meshTriangleIndices);

    for (Mesh &mesh : container->meshes)
    {
        mesh.SetPositionsHandle(
            {internStates.positions.ResolveIndex(mesh.GetPositionsHandle().index)});
        mesh.SetNormalsHandle({internStates.normals.ResolveIndex(mesh.GetNormalsHandle().index)});
        mesh.SetUVsHandle({internStates.uvs.ResolveIndex(mesh.GetUVsHandle().index)});
        mesh.SetUVIndicesHandle(
            {internStates.uvIndices.ResolveIndex(mesh.GetUVIndicesHandle().index)});
        mesh.SetIndicesHandle({internStates.indices.ResolveIndex(mesh.GetIndicesHandle().index)});
    }
#endif

    // Get all materials
    for (size_t sceneIndex = 0; sceneIndex < traversalDatas.size(); sceneIndex++)
    {
        USDPrimTraversalData &data = traversalDatas[sceneIndex];
        size_t meshOffset = sceneMeshOffsets[sceneIndex];

        for (size_t i = 0; i < data.meshes.size(); i++)
        {
            const size_t meshIndex = meshOffset + i;
            Mesh &mesh = container->meshes[meshIndex];
            pxr::UsdPrim prim = data.meshes[i];

            pxr::UsdShadeMaterialBindingAPI bindingApi(prim);
            std::vector<pxr::UsdGeomSubset> subsets = bindingApi.GetMaterialBindSubsets();

            if (subsets.size() == 0)
            {
                pxr::UsdShadeMaterial material = GetPrimMaterial(prim);

                if (!material)
                {
                    printf("%s is missing a material\n", prim.GetPath().GetString().c_str());
                    continue;
                }

                mesh.SetMaterialIndex(AddMaterialToMap(materialMap, materials, material));
            }
            else
            {
                Array<int> faceMaterialIndices;

                pxr::VtIntArray faceCounts;
                pxr::UsdGeomMesh usdMesh(prim);
                bool result = usdMesh.GetFaceVertexCountsAttr().Get(&faceCounts, timeCode);
                YBI_ASSERT(result);

                GetFaceMaterialIndices(
                    faceMaterialIndices, faceCounts, subsets, materialMap, materials, timeCode);
#ifdef USD_CONCURRENT_MESH_LOAD
                size_t provisionalIndex = internStates.faceMaterialIndices.Intern(
                    scratch.arena, std::move(faceMaterialIndices), meshIndex);
                mesh.SetFaceMaterialIndicesHandle({static_cast<uint32_t>(provisionalIndex)});
#else
                BufferHandle faceMaterialHandle{
                    static_cast<uint32_t>(container->meshFaceMaterialIndices.size())};
                container->meshFaceMaterialIndices.push_back(std::move(faceMaterialIndices));
                mesh.SetFaceMaterialIndicesHandle(faceMaterialHandle);
#endif
            }
        }
    }

#ifdef USD_CONCURRENT_MESH_LOAD
    internStates.faceMaterialIndices.MoveSortedBuffersTo(&container->meshFaceMaterialIndices);
    for (Mesh &mesh : container->meshes)
    {
        mesh.SetFaceMaterialIndicesHandle({internStates.faceMaterialIndices.ResolveIndex(
            mesh.GetFaceMaterialIndicesHandle().index)});
    }
#endif

    // Handle materials
    bool dumpedDebugMaterial = false;
    const char *debugMaterialFilter = std::getenv("YBI_DEBUG_MATERIAL");

    for (pxr::UsdShadeMaterial &material : materials)
    {
        MaterialHandle &handle = container->materialHandles.emplace_back();

        pxr::UsdShadeShader shader = GetMaterialSurfaceShader(material);

        if (shader)
        {
            const std::string materialPath = material.GetPath().GetString();

            pxr::TfToken token;
            if (shader.GetShaderId(&token) && token == pxr::TfToken("UsdPreviewSurface"))
            {
                handle.isGraph = false;
                handle.index = container->materials.size();

                Material &mat = container->materials.emplace_back();

                mat.path = material.GetPath().GetString();
                for (const pxr::UsdShadeInput &input : shader.GetInputs())
                {
                    if (input.GetBaseName().GetString() != "diffuseColor")
                        continue;

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

                    ParseSTInput(shader);

                    pxr::SdfAssetPath assetPath;
                    bool result = fileInput.Get(&assetPath);
                    YBI_ASSERT(result);

                    ParseUdimPaths(assetPath, fileInput.GetAttr(), timeCode, mat.diffuseUdimPaths);

#if 0
                    pxr::TfToken wrapS, wrapT;
                    pxr::UsdShadeInput wrapSInput = shader.GetInput(pxr::TfToken("wrapS"));
                    pxr::UsdShadeInput wrapTInput = shader.GetInput(pxr::TfToken("wrapT"));

                    if (wrapSInput && wrapSInput.Get(&wrapS))
                    {
                        printf("wraps: %s\n", wrapS.GetText());
                    }
                    if (wrapTInput && wrapTInput.Get(&wrapT))
                    {
                        printf("wrapt: %s\n", wrapT.GetText());
                    }
#endif
                }
            }
            else if (token == pxr::TfToken("ND_open_pbr_surface_surfaceshader"))
            {
                handle.isGraph = true;
                handle.index = container->openPbrMaterials.size();

                OpenPbrMaterial &mat = container->openPbrMaterials.emplace_back();
                LoadOpenPbrMaterial(shader, timeCode, &mat);
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
            // YBI_ASSERT(0);
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

    // Handle lights
    for (size_t sceneIndex = 0; sceneIndex < traversalDatas.size(); sceneIndex++)
    {
        const USDPrimTraversalData &data = traversalDatas[sceneIndex];

        for (size_t i = 0; i < data.rectLights.size(); i++)
        {
            pxr::UsdLuxRectLight rectLight(data.rectLights[i]);

            pxr::UsdAttribute widthAttr = rectLight.GetWidthAttr();
            pxr::UsdAttribute heightAttr = rectLight.GetHeightAttr();
            pxr::UsdAttribute intensityAttr = rectLight.GetIntensityAttr();
            pxr::UsdAttribute exposureAttr = rectLight.GetExposureAttr();
            pxr::UsdAttribute colorAttr = rectLight.GetColorAttr();
            pxr::UsdAttribute normalizeAttr = rectLight.GetNormalizeAttr();

            float width, height, intensity;
            float exposure;
            bool normalize;
            pxr::GfVec3f radiance;
            widthAttr.Get(&width, timeCode);
            heightAttr.Get(&height, timeCode);
            intensityAttr.Get(&intensity, timeCode);
            exposureAttr.Get(&exposure, timeCode);
            colorAttr.Get(&radiance, timeCode);
            normalizeAttr.Get(&normalize, timeCode);

            const pxr::GfMatrix4d localToWorld =
                xformCache.GetLocalToWorldTransform(rectLight.GetPrim());

            AffineTransform transform = ConvertAffineTransform(localToWorld);

            printf("chungus: %f %f %f %f\n", width, height, intensity, exposure);

            RectLight light;
            light.localToWorld = transform;
            light.width = width;
            light.height = height;
            light.radiance =
                Vec3f(radiance[0], radiance[1], radiance[2]) * intensity * exp2f(exposure);

            if (normalize)
            {
                const Vec3f worldWidth =
                    transform.TransformVector(Vec3f(width, 0.f, 0.f));
                const Vec3f worldHeight =
                    transform.TransformVector(Vec3f(0.f, height, 0.f));
                const float worldArea = Length(Cross(worldWidth, worldHeight));

                if (worldArea > 0.f)
                {
                    light.radiance /= worldArea;
                }
                else
                {
                    light.radiance = Vec3f(0.f);
                }
            }

            container->rectLights.push_back(light);

            // pxr::UsdAttribute intensityAttr = rectLight.GetIntensityAttr();
            // intensityAttr.GetSamp
        }
    }

    // Defer loading of geometry attributes until all material graphs have been traversed.
    {
        for (size_t sceneIndex = 0; sceneIndex < traversalDatas.size(); sceneIndex++)
        {
            USDPrimTraversalData &data = traversalDatas[sceneIndex];
            size_t meshOffset = sceneMeshOffsets[sceneIndex];
            for (size_t i = 0; i < data.meshes.size(); i++)
            {
                const size_t meshIndex = meshOffset + i;
                Mesh &mesh = container->meshes[meshIndex];
                pxr::UsdPrim prim = data.meshes[i];

                // TODO: should make it so that only graphs can enter this path too?
                if (mesh.HasMaterial())
                {
                    const pxr::SdfPath primPath = prim.GetPath();
                    auto loadAttribute = [stage, primPath, timeCode](const std::string &propName,
                                                                     MeshAttribute &outAttribute) {
                        return LoadPrimvar(stage, primPath, timeCode, propName, outAttribute);
                    };

                    // printf("size of lambda: %zu\n", sizeof(loadAttribute));

                    importContext->pendingMeshes.push_back({meshIndex, std::move(loadAttribute)});
                }
            }
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
                    YBI_LOGFATAL("Geometry type: %s not supported\n",
                                 ConvertGeometryTypeToString(instance.handle.type));
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
