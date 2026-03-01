#include "tessellation/subdivision.h"

#include "tessellation/edge_map_validation.h"
#include "tessellation/subdivision_patch_types.h"
#include "util/assert.h"

#include <opensubdiv/far/patchMap.h>
#include <opensubdiv/far/patchTableFactory.h>
#include <opensubdiv/far/primvarRefiner.h>
#include <opensubdiv/far/ptexIndices.h>
#include <opensubdiv/far/topologyDescriptor.h>
#include <opensubdiv/far/topologyRefinerFactory.h>
#include <pxr/base/gf/vec2f.h>

#include <algorithm>
#include <cmath>
#include <cstring>
#include <cstdint>
#include <cstdio>
#include <limits>
#include <string>
#include <unordered_map>
#include <vector>

using namespace OpenSubdiv;

YBI_NAMESPACE_BEGIN

template <typename T> struct LimitEvalValueTraits;

template <> struct LimitEvalValueTraits<float>
{
    static float Zero()
    {
        return 0.0f;
    }
};

template <> struct LimitEvalValueTraits<float2>
{
    static float2 Zero()
    {
        return make_float2(0.0f);
    }
};

template <> struct LimitEvalValueTraits<float3>
{
    static float3 Zero()
    {
        return make_float3(0.0f);
    }
};

template <> struct LimitEvalValueTraits<float4>
{
    static float4 Zero()
    {
        return make_float4(0.0f, 0.0f, 0.0f, 0.0f);
    }
};

template <typename T> struct LimitEvalValue
{
    T value = LimitEvalValueTraits<T>::Zero();

    void Clear()
    {
        value = LimitEvalValueTraits<T>::Zero();
    }

    void AddWithWeight(const LimitEvalValue &src, float w)
    {
        value += src.value * w;
    }
};

using LimitEvalFloat = LimitEvalValue<float>;
using LimitEvalFloat2 = LimitEvalValue<float2>;
using LimitEvalFloat3 = LimitEvalValue<float3>;
using LimitEvalFloat4 = LimitEvalValue<float4>;
using LimitEvalVertex = LimitEvalFloat3;
using LimitEvalFVar2 = LimitEvalFloat2;

struct SubdivisionLimitEvalAttribute
{
    std::string name;
    AttributeType type = AttributeType::Unknown;
    PrimvarInterpolation interpolation = PrimvarInterpolation::Unknown;
    int fvarChannel = -1;

    std::vector<LimitEvalFloat> valuesFloat;
    std::vector<LimitEvalFloat2> valuesFloat2;
    std::vector<LimitEvalFloat3> valuesFloat3;
    std::vector<LimitEvalFloat4> valuesFloat4;
};

template <typename T> static Array<T> BytesToArray(const Array<uint8_t> &bytes)
{
    YBI_ASSERT((bytes.size() % sizeof(T)) == 0);
    const size_t count = bytes.size() / sizeof(T);
    Array<T> out(count);
    if (count > 0)
    {
        memcpy(out.data(), bytes.data(), bytes.size());
    }
    return out;
}

// clang-format off
#include "tessellation/tessellation_adaptive_limit_eval.h"
#include "tessellation/tessellation_adaptive_edge_ops.h"
#include "tessellation/tessellation_adaptive_patch_build.h"
#include "tessellation/tessellation_adaptive_util.h"
#include "tessellation/tessellation_adaptive_patch_quad_write.h"
#include "tessellation/tessellation_adaptive_mesh_build.h"
// clang-format on

bool SubdivideAdaptive(const SubdivisionMesh &mesh,
                      const SubdivisionRunOptions &options,
                      SubdivisionRunResult *outResult)
{
    if (!outResult)
    {
        return false;
    }

    CreasePairs creases = BuildCreasePairs(mesh);
    Far::TopologyDescriptor d = {};
    d.numVertices = int(mesh.vertices.size());
    d.numFaces = int(mesh.vertsPerFace.size());
    d.numVertsPerFace = mesh.vertsPerFace.data();
    d.vertIndicesPerFace = mesh.indices.data();
    d.numCreases = int(creases.weights.size());
    d.creaseVertexIndexPairs = creases.pairs.data();
    d.creaseWeights = creases.weights.data();
    d.numCorners = int(mesh.cornerIndices.size());
    d.cornerVertexIndices = mesh.cornerIndices.data();
    d.cornerWeights = mesh.cornerSharpnesses.data();
    d.numHoles = int(mesh.holeIndices.size());
    d.holeIndices = mesh.holeIndices.data();
    std::vector<const Attribute *> tessellationAttributes;
    tessellationAttributes.reserve(mesh.attributes.size());
    for (const Attribute &attr : mesh.attributes)
    {
        if (!AttributeTypeIsFloat(attr.type))
        {
            continue;
        }
        if (attr.interpolation == PrimvarInterpolation::Vertex ||
            attr.interpolation == PrimvarInterpolation::Varying ||
            attr.interpolation == PrimvarInterpolation::FaceVarying)
        {
            tessellationAttributes.push_back(&attr);
        }
    }

    std::vector<Far::TopologyDescriptor::FVarChannel> fvarChannels;
    std::vector<Array<int>> generatedFVarIndices;
    std::unordered_map<const Attribute *, int> fvarChannelByAttribute;
    for (const Attribute *attr : tessellationAttributes)
    {
        if (!attr || attr->interpolation != PrimvarInterpolation::FaceVarying)
        {
            continue;
        }
        const size_t elementSize = AttributeTypeGetSize(attr->type);
        if (elementSize == 0 || (attr->data.size() % elementSize) != 0)
        {
            continue;
        }
        const int numValues = int(attr->data.size() / elementSize);
        if (numValues <= 0)
        {
            continue;
        }

        const int *valueIndices = nullptr;
        if (attr->indices.size() == mesh.indices.size())
        {
            valueIndices = attr->indices.data();
        }
        else if (size_t(numValues) == mesh.indices.size())
        {
            generatedFVarIndices.emplace_back(mesh.indices.size());
            for (size_t i = 0; i < mesh.indices.size(); ++i)
            {
                generatedFVarIndices.back()[i] = int(i);
            }
            valueIndices = generatedFVarIndices.back().data();
        }
        else
        {
            continue;
        }
        fvarChannelByAttribute[attr] = int(fvarChannels.size());
        Far::TopologyDescriptor::FVarChannel channel = {};
        channel.numValues = numValues;
        channel.valueIndices = valueIndices;
        fvarChannels.push_back(channel);
    }
    if (!fvarChannels.empty())
    {
        d.numFVarChannels = int(fvarChannels.size());
        d.fvarChannels = fvarChannels.data();
    }

    Sdc::Options sdcOptions;
    switch (mesh.interpolationRule)
    {
        case BOUNDARY_INTERPOLATION_NONE:
            sdcOptions.SetVtxBoundaryInterpolation(Sdc::Options::VTX_BOUNDARY_NONE);
            break;
        case BOUNDARY_INTERPOLATION_EDGE:
            sdcOptions.SetVtxBoundaryInterpolation(Sdc::Options::VTX_BOUNDARY_EDGE_ONLY);
            break;
        case BOUNDARY_INTERPOLATION_EDGE_AND_CORNER:
        default:
            sdcOptions.SetVtxBoundaryInterpolation(Sdc::Options::VTX_BOUNDARY_EDGE_AND_CORNER);
            break;
    }
    switch (mesh.fvarLinearInterpolation)
    {
        case FVAR_LINEAR_NONE:
            sdcOptions.SetFVarLinearInterpolation(Sdc::Options::FVAR_LINEAR_NONE);
            break;
        case FVAR_LINEAR_CORNERS_ONLY:
            sdcOptions.SetFVarLinearInterpolation(Sdc::Options::FVAR_LINEAR_CORNERS_ONLY);
            break;
        case FVAR_LINEAR_CORNERS_PLUS1:
            sdcOptions.SetFVarLinearInterpolation(Sdc::Options::FVAR_LINEAR_CORNERS_PLUS1);
            break;
        case FVAR_LINEAR_CORNERS_PLUS2:
            sdcOptions.SetFVarLinearInterpolation(Sdc::Options::FVAR_LINEAR_CORNERS_PLUS2);
            break;
        case FVAR_LINEAR_BOUNDARIES:
            sdcOptions.SetFVarLinearInterpolation(Sdc::Options::FVAR_LINEAR_BOUNDARIES);
            break;
        case FVAR_LINEAR_ALL:
            sdcOptions.SetFVarLinearInterpolation(Sdc::Options::FVAR_LINEAR_ALL);
            break;
        default:
            sdcOptions.SetFVarLinearInterpolation(Sdc::Options::FVAR_LINEAR_CORNERS_PLUS1);
            break;
    }
    sdcOptions.SetCreasingMethod(ToSdcCreasingMethod(options.creasingMethod));
    sdcOptions.SetTriangleSubdivision(ToSdcTriangleSubdivision(mesh.triangleSubdivisionRule));

    Far::TopologyRefinerFactory<Far::TopologyDescriptor>::Options o(
        SchemeFromString(options.subdivisionScheme), sdcOptions);

    SubdivisionEdgeMap edgeMap;
    edgeMap.reserve(mesh.indices.size());
    const int edgesWithOver2Faces = CountNonManifoldEdges(edgeMap);

    Far::TopologyRefiner *refiner =
        Far::TopologyRefinerFactory<Far::TopologyDescriptor>::Create(d, o);
    if (!refiner)
    {
        return false;
    }
    const bool hasRefinerFVar = refiner->GetNumFVarChannels() > 0;
    YBI_ERROR(fvarChannels.empty() || hasRefinerFVar,
              "Subdivision face-varying attrs present but refiner has no fvar channels\n");

    int nextGeneratedVertexId = int(mesh.vertices.size());
    const std::vector<SubdivisionPatch> patches =
        BuildSubdivisionPatches(mesh, *refiner, edgeMap, nextGeneratedVertexId);

    const EdgeMapChecks edgeChecks = RunEdgeMapChecks(mesh, patches, edgeMap);
    if (!edgeChecks.ok)
    {
        delete refiner;
        return false;
    }

    const bool enableFVarTables = !fvarChannels.empty() && hasRefinerFVar;

    Far::TopologyRefiner::AdaptiveOptions adaptiveOptions(options.level);
    adaptiveOptions.considerFVarChannels = enableFVarTables;
    refiner->RefineAdaptive(adaptiveOptions);

    Far::PatchTableFactory::Options patchOptions(options.level);
    patchOptions.endCapType = Far::PatchTableFactory::Options::ENDCAP_GREGORY_BASIS;
    patchOptions.useInfSharpPatch = (d.numCreases > 0) || (d.numCorners > 0);
    patchOptions.generateFVarTables = enableFVarTables;
    patchOptions.generateFVarLegacyLinearPatches = false;
    const Far::PatchTable *patchTable = Far::PatchTableFactory::Create(*refiner, patchOptions);
    if (!patchTable)
    {
        delete refiner;
        return false;
    }

    Far::PatchMap patchMap(*patchTable);
    const std::vector<LimitEvalVertex> limitValues =
        BuildLimitEvalVertices(*refiner, *patchTable, mesh.vertices);
    std::vector<SubdivisionLimitEvalAttribute> limitEvalAttributes;
    limitEvalAttributes.reserve(tessellationAttributes.size());
    for (const Attribute *attr : tessellationAttributes)
    {
        if (!attr)
        {
            continue;
        }
        SubdivisionLimitEvalAttribute outAttr = {};
        outAttr.name = attr->name;
        outAttr.type = attr->type;
        outAttr.interpolation = attr->interpolation;
        if (attr->interpolation == PrimvarInterpolation::FaceVarying)
        {
            auto found = fvarChannelByAttribute.find(attr);
            if (found == fvarChannelByAttribute.end() || found->second < 0)
            {
                continue;
            }
            outAttr.fvarChannel = found->second;
        }

        bool ok = false;
        if (attr->type == AttributeType::Float)
        {
            Array<float> coarse = BytesToArray<float>(attr->data);
            outAttr.valuesFloat = BuildLimitEvalValues<LimitEvalFloat, float>(
                *refiner, *patchTable, coarse, attr->interpolation, outAttr.fvarChannel);
            ok = !outAttr.valuesFloat.empty();
        }
        else if (attr->type == AttributeType::Float2)
        {
            Array<float2> coarse = BytesToArray<float2>(attr->data);
            outAttr.valuesFloat2 = BuildLimitEvalValues<LimitEvalFloat2, float2>(
                *refiner, *patchTable, coarse, attr->interpolation, outAttr.fvarChannel);
            ok = !outAttr.valuesFloat2.empty();
        }
        else if (attr->type == AttributeType::Float3)
        {
            Array<float3> coarse = BytesToArray<float3>(attr->data);
            outAttr.valuesFloat3 = BuildLimitEvalValues<LimitEvalFloat3, float3>(
                *refiner, *patchTable, coarse, attr->interpolation, outAttr.fvarChannel);
            ok = !outAttr.valuesFloat3.empty();
        }
        else if (attr->type == AttributeType::Float4)
        {
            Array<float4> coarse = BytesToArray<float4>(attr->data);
            outAttr.valuesFloat4 = BuildLimitEvalValues<LimitEvalFloat4, float4>(
                *refiner, *patchTable, coarse, attr->interpolation, outAttr.fvarChannel);
            ok = !outAttr.valuesFloat4.empty();
        }
        if (ok)
        {
            limitEvalAttributes.push_back(std::move(outAttr));
        }
    }

    const float3 eye = options.eye;
    const float3 lookAt = options.lookAt;

    int tmaxComputedEdges = 0;
    const std::vector<SubdivisionPatch> splitPatches =
        DiagSplitPatches(patches,
                         edgeMap,
                         nextGeneratedVertexId,
                         options.maxDiagSplitDepth,
                         patchMap,
                         *patchTable,
                         limitValues,
                         options.sampleSteps,
                         options.pixelSpacing,
                         options.splitThreshold,
                         eye,
                         lookAt,
                         options.viewportWidth,
                         options.viewportHeight,
                         options.verticalFovDegrees,
                         options.useCameraMatrices,
                         options.cameraFromWorld,
                         options.clipFromCamera,
                         &tmaxComputedEdges);

    int missingStoredPatchParams = 0;
    int badUniformFactor = 0;
    if (!VerifyInitializedUniformEdgesHaveStoredPatchParams(
            edgeMap, &missingStoredPatchParams, &badUniformFactor))
    {
        delete patchTable;
        delete refiner;
        return false;
    }

    int patchQuadVerts = 0;
    int patchQuadCount = 0;
    if (!options.patchQuadObjPath.empty())
    {
        if (!WriteLeafPatchCornerQuadsObj(splitPatches,
                                          edgeMap,
                                          patchMap,
                                          *patchTable,
                                          limitValues,
                                          options.patchQuadObjPath,
                                          &patchQuadVerts,
                                          &patchQuadCount))
        {
            delete patchTable;
            delete refiner;
            return false;
        }
    }

    int innerGridVerts = 0;
    int innerGridTris = 0;
    int innerGridOnlyTris = 0;
    int stitchingOnlyTris = 0;
    std::vector<int> *triPatchFaceIds = nullptr;
    std::vector<int> *triCoarseFaceIds = nullptr;
    std::vector<int> *triPtexFaceIds = nullptr;
    std::vector<int> *triQuadrants = nullptr;
    if (options.generateTriangleMetadata)
    {
        triPatchFaceIds = &outResult->trianglePatchFaceIds;
        triCoarseFaceIds = &outResult->triangleCoarseFaceIds;
        triPtexFaceIds = &outResult->trianglePtexFaceIds;
        triQuadrants = &outResult->triangleQuadrants;
    }
    else
    {
        outResult->trianglePatchFaceIds.clear();
        outResult->triangleCoarseFaceIds.clear();
        outResult->trianglePtexFaceIds.clear();
        outResult->triangleQuadrants.clear();
    }
    if (!BuildLeafPatchStitchedMesh(splitPatches,
                                    edgeMap,
                                    patchMap,
                                    *patchTable,
                                    limitValues,
                                    limitEvalAttributes,
                                    nextGeneratedVertexId,
                                    &outResult->mesh,
                                    triPatchFaceIds,
                                    triCoarseFaceIds,
                                    triPtexFaceIds,
                                    triQuadrants,
                                    &innerGridVerts,
                                    &innerGridTris,
                                    &innerGridOnlyTris,
                                    &stitchingOnlyTris))
    {
        delete patchTable;
        delete refiner;
        return false;
    }

    outResult->refinedMaxLevel = refiner->GetMaxLevel();
    outResult->totalPatches = patchTable->GetNumPatchesTotal();
    outResult->subdivisionPatchCount = patches.size();
    outResult->diagSplitPatchCount = splitPatches.size();
    outResult->generatedVertexCount = nextGeneratedVertexId - int(mesh.vertices.size());
    outResult->midpointEdges = CountEdgesWithMidpointVertex(edgeMap);
    outResult->edgeTMaxComputed = tmaxComputedEdges;
    outResult->totalComputedEdges = CountEdgesWithComputedTMax(edgeMap);
    outResult->patchQuadVerts = patchQuadVerts;
    outResult->patchQuadCount = patchQuadCount;
    outResult->innerGridTriangleCount = innerGridOnlyTris;
    outResult->stitchingTriangleCount = stitchingOnlyTris;
    outResult->controlCageUniqueEdges = edgeMap.size();
    outResult->boundaryEdges = CountBoundaryEdges(edgeMap);
    outResult->controlCageEdgesWithOver2Faces = edgesWithOver2Faces;

    delete patchTable;
    delete refiner;
    return true;
}

YBI_NAMESPACE_END
