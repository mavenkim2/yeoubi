#include "scene/subdivision.h"
#include "scene/attributes.h"
#include "scene/scene.h"
#include "util/assert.h"
#include "util/base.h"
#include "util/float2.h"
#include "util/float4x4.h"
#include "util/host_memory_arena.h"
#include <algorithm>
#include <opensubdiv/far/patchMap.h>
#include <opensubdiv/far/patchTable.h>
#include <opensubdiv/far/patchTableFactory.h>
#include <opensubdiv/far/primvarRefiner.h>
#include <opensubdiv/far/ptexIndices.h>
#include <opensubdiv/far/topologyDescriptor.h>
#include <opensubdiv/sdc/types.h>
#include <opensubdiv/version.h>
#include <vector>

YBI_NAMESPACE_BEGIN

namespace Far = OpenSubdiv::Far;
namespace Sdc = OpenSubdiv::Sdc;

template <typename T>
struct OsdData
{
    T value;

    void Clear()
    {
        memset(&value, 0, sizeof(T));
    }
    void AddWithWeight(const OsdData &src, float w)
    {
        value += src.value * w;
    }
};

static OpenSubdiv::Sdc::Options::VtxBoundaryInterpolation ToOsdVtxBoundary(BoundaryInterpolation r)
{
    switch (r)
    {
        case BOUNDARY_INTERPOLATION_NONE:
            return Sdc::Options::VTX_BOUNDARY_NONE;
        case BOUNDARY_INTERPOLATION_EDGE:
            return Sdc::Options::VTX_BOUNDARY_EDGE_ONLY;
        case BOUNDARY_INTERPOLATION_EDGE_AND_CORNER:
            return Sdc::Options::VTX_BOUNDARY_EDGE_AND_CORNER;
        default:
            return Sdc::Options::VTX_BOUNDARY_EDGE_AND_CORNER;
    }
}

static OpenSubdiv::Sdc::Options::FVarLinearInterpolation ToOsdFVarLinear(FVarLinearInterpolation f)
{
    switch (f)
    {
        case FVAR_LINEAR_NONE:
            return Sdc::Options::FVAR_LINEAR_NONE;
        case FVAR_LINEAR_CORNERS_ONLY:
            return Sdc::Options::FVAR_LINEAR_CORNERS_ONLY;
        case FVAR_LINEAR_CORNERS_PLUS1:
            return Sdc::Options::FVAR_LINEAR_CORNERS_PLUS1;
        case FVAR_LINEAR_CORNERS_PLUS2:
            return Sdc::Options::FVAR_LINEAR_CORNERS_PLUS2;
        case FVAR_LINEAR_BOUNDARIES:
            return Sdc::Options::FVAR_LINEAR_BOUNDARIES;
        case FVAR_LINEAR_ALL:
            return Sdc::Options::FVAR_LINEAR_ALL;
        default:
            return Sdc::Options::FVAR_LINEAR_CORNERS_ONLY;
    }
}

template <typename DataType, typename T>
static MemoryView<OsdData<DataType>> InterpolateVertex(HostMemoryArena &arena,
                                                       const T &array,
                                                       const Far::TopologyRefiner *refiner,
                                                       const Far::PatchTable *patchTable)
{

    const int numLevels = refiner->GetNumLevels();
    const int numVertices = refiner->GetNumVerticesTotal();
    const int numLocalPoints = patchTable->GetNumLocalPoints();

    YBI_ASSERT(array.size() == refiner->GetLevel(0).GetNumVertices());

    size_t typeSize = sizeof(DataType);
    size_t numFloats = typeSize / sizeof(float);

    MemoryView<OsdData<DataType>> view =
        arena.PushArray<OsdData<DataType>>(numVertices + numLocalPoints);

    memcpy(view.data(), array.data(), array.size() * typeSize);

    OsdData<DataType> *src = view.data();

    for (int level = 1; level < numLevels; level++)
    {
        OsdData<DataType> *dst = src + refiner->GetLevel(level - 1).GetNumVertices();
        Far::PrimvarRefiner(*refiner).Interpolate(level, src, dst);
        src = dst;
    }
    patchTable->ComputeLocalPointValues(view.data(), view.data() + numVertices);
    return view;
}

#if 0
template <typename T>
static void InterpolateVarying(HostMemoryArena &arena,
                               const T &array,
                               const Far::TopologyRefiner *refiner,
                               const Far::PatchTable *patchTable,
                               int channel)
{
    const int numVerticesTotal = refiner->GetNumFVarValuesTotal(channel);
    const int numLocalPointsTotal = patchTable->GetNumLocalPointsVarying();

    MemoryView<OsdData> view = arena.PushArray<OsdData>(numVerticesTotal + numLocalPointsTotal);

    const int numLevels = refiner->GetNumLevels();
    OsdData *src = view.data();

    for (int level = 1; level < numLevels; level++)
    {
        OsdData *dst = src + refiner->GetLevel(level - 1).GetNumFVarValues(channel);
        Far::PrimvarRefiner(*refiner).InterpolateVarying(level, src, dst);
        src = dst;
    }
    patchTable->ComputeLocalPointValuesVarying(src, view.data() + numVerticesTotal);
}

template <typename T>
static void InterpolateFaceVarying(MemoryView<OsdData> &view,
                                   const T &array,
                                   const Far::TopologyRefiner *refiner,
                                   const Far::PatchTable *patchTable,
                                   int channel)
{
    const int numVerticesTotal = refiner->GetNumFVarValuesTotal(channel);
    const int numLocalPointsTotal = patchTable->GetNumLocalPointsFaceVarying();
    const int numLevels = refiner->GetNumLevels();
    OsdData *src = view.data();

    for (int level = 1; level < numLevels; level++)
    {
        OsdData *dst = src + refiner->GetLevel(level - 1).GetNumFVarValues(channel);
        Far::PrimvarRefiner(*refiner).InterpolateFaceVarying(level, src, dst);
        src = dst;
    }
    patchTable->ComputeLocalPointValuesVarying(src, view.data() + numVerticesTotal);
}
#endif

struct LimitSurfaceSample
{
    float2 uv;
    int faceID;
};

static void GenerateLimitSurfaceSamples(const Far::PatchMap *patchMap,
                                        const Far::PatchTable *patchTable)
{
    // for each face
    // get the edge rates of each edge
    // this determines the internal grid
    // determine the edges (if any) that have different edge rates from their opposites

    for (int patch = 0; patch < patchTable->GetNumPatchesTotal(); patch++)
    {
        // const Far::PatchTable::PatchHandle *handle = patchMap->FindPatch(fid);
        // patchTable->GetPatchVertices()
    }
    // patchTable->GetPatchVert
}

// https://techmatt.github.io/pdfs/diagSplit.pdf
namespace DiagSplit
{

static const int NON_UNIFORM = -1;
static const float2 CORNER_TO_UV[4] = {
    make_float2(0, 0),
    make_float2(1, 0),
    make_float2(1, 1),
    make_float2(0, 1),
};
static const float2 UV_DELTA[4] = {
    make_float2(1, 0),
    make_float2(0, 1),
    make_float2(-1, 0),
    make_float2(0, -1),
};

struct DiagSplitParams
{
    const Far::PatchMap *patchMap;
    const Far::PatchTable *patchTable;
    const MemoryView<OsdData<float3>> positions;
    int N;   // number of times to sample edge in T()
    float R; // vertex pixel spacing goal (in pixels)
    int splitThreshold;
    float4x4 viewProj; // row-major world/camera -> clip. ToScreen uses this to project to pixels.
    int viewportWidth;
    int viewportHeight;
};

struct SubPatch
{
    float2 parametricCorners[4];
    int edgeRates[4];
    int ptexFace;

    SubPatch(int face) : ptexFace(face)
    {
        edgeRates[0] = DiagSplit::NON_UNIFORM;
        edgeRates[1] = DiagSplit::NON_UNIFORM;
        edgeRates[2] = DiagSplit::NON_UNIFORM;
        edgeRates[3] = DiagSplit::NON_UNIFORM;

        parametricCorners[0] = make_float2(0.f, 0.f);
        parametricCorners[1] = make_float2(1.f, 0.f);
        parametricCorners[2] = make_float2(1.f, 1.f);
        parametricCorners[3] = make_float2(0.f, 1.f);
    }
};

__forceinline static int Next(int edge, int numFaceVertices)
{
    return (edge + 1) % numFaceVertices;
}

static float3 EvaluatePosition(const DiagSplitParams &params, int fid, float2 uv)
{
    const Far::PatchTable::PatchHandle *handle = params.patchMap->FindPatch(fid, uv.x, uv.y);
    auto cvIndices = params.patchTable->GetPatchVertices(*handle);

    float pWeights[20];
    params.patchTable->EvaluateBasis(*handle, uv.x, uv.y, pWeights);
    float3 p = make_float3(0.f);
    for (int cv = 0; cv < cvIndices.size(); cv++)
    {
        p += params.positions[cvIndices[cv]].value * pWeights[cv];
    }
    return p;
}

static float2 ToScreen(const DiagSplitParams &params, const float3 &p)
{
    float4 clip = mul(params.viewProj, make_float4(p.x, p.y, p.z, 1.f));
    if (clip.w <= 0.f)
        return make_float2(0.f);
    float invW = 1.f / clip.w;
    float ndcX = clip.x * invW;
    float ndcY = clip.y * invW;
    float sx = (ndcX * 0.5f + 0.5f) * (float)params.viewportWidth;
    float sy = (1.f - (ndcY * 0.5f + 0.5f)) * (float)params.viewportHeight;
    return make_float2(sx, sy);
}

static int T(const DiagSplitParams &params, const float2 &uvStart, const float2 &uvEnd, int fid)
{
    const int N = params.N;
    const float R = params.R;

    YBI_ASSERT(N > 1);

    float maxLi = 0.f;
    float sumLi = 0.f;

    float3 pPrev = EvaluatePosition(params, fid, uvStart);
    float2 screenPrev = ToScreen(params, pPrev);

    for (int i = 1; i < N; i++)
    {
        float2 uv = lerp(uvStart, uvEnd, float(i) / (N - 1));
        float3 p = EvaluatePosition(params, fid, uv);
        float2 screenP = ToScreen(params, p);
        float Li = length(screenP - screenPrev);
        sumLi += Li;
        maxLi = std::max(maxLi, Li);
        pPrev = p;
        screenPrev = screenP;
    }

    int tMin = (int)ceilf(sumLi / R);
    int tMax = (int)ceilf((N - 1) * maxLi / R);

    if (tMax - tMin >= params.splitThreshold)
    {
        return DiagSplit::NON_UNIFORM;
    }
    return tMax;
}

static float2 PartitionEdge(const DiagSplitParams &params,
                            int fid,
                            float2 uvStart,
                            float2 uvEnd,
                            int edgeFactor,
                            int &t0,
                            int &t1)
{
    if (edgeFactor == DiagSplit::NON_UNIFORM)
    {
        float2 uv = (uvStart + uvEnd) / 2.f;

        t0 = T(params, uvStart, uv, fid);
        t1 = T(params, uv, uvEnd, fid);
        return uv;
    }
    else
    {
        t0 = (int)floorf(edgeFactor / 2.f);
        t1 = edgeFactor - t0;

        float2 uv = lerp(uvStart, uvEnd, float(t0) / edgeFactor);
        return uv;
    }
}

static void Split(const DiagSplitParams &params, SubPatch &patch)
{
    const int numFaceVertices = 4;
    for (int edge = 0; edge < numFaceVertices; edge++)
    {
        if (patch.edgeRates[edge] == DiagSplit::NON_UNIFORM)
        {
            int next = Next(edge, numFaceVertices);
            int opp = Next(edge + 1, numFaceVertices);
            int prev = Next(edge + 2, numFaceVertices);

            float2 uvStart = patch.parametricCorners[edge];
            float2 uvEnd = patch.parametricCorners[next];
            int edgeFactor = patch.edgeRates[edge];

            int edgeT0, edgeT1;
            float2 uvMid =
                PartitionEdge(params, patch.ptexFace, uvStart, uvEnd, edgeFactor, edgeT0, edgeT1);

            uvStart = patch.parametricCorners[opp];
            uvEnd = patch.parametricCorners[prev];
            edgeFactor = patch.edgeRates[opp];
            int edgeOppT0, edgeOppT1;
            float2 uvMidOpp = PartitionEdge(
                params, patch.ptexFace, uvStart, uvEnd, edgeFactor, edgeOppT0, edgeOppT1);

            SubPatch newPatch0(patch.ptexFace);
            newPatch0.parametricCorners[0] = patch.parametricCorners[edge];
            newPatch0.parametricCorners[1] = uvMid;
            newPatch0.parametricCorners[2] = uvMidOpp;
            newPatch0.parametricCorners[3] = patch.parametricCorners[prev];
            newPatch0.edgeRates[0] = edgeT0;
            newPatch0.edgeRates[1] = patch.edgeRates[next];
            newPatch0.edgeRates[2] = edgeOppT1;
            newPatch0.edgeRates[3] = patch.edgeRates[prev];

            SubPatch newPatch1(patch.ptexFace);
            newPatch1.ptexFace = patch.ptexFace;
            newPatch1.parametricCorners[0] = uvMid;
            newPatch1.parametricCorners[1] = patch.parametricCorners[next];
            newPatch1.parametricCorners[2] = patch.parametricCorners[opp];
            newPatch1.parametricCorners[3] = uvMidOpp;
            newPatch1.edgeRates[0] = edgeT1;
            newPatch1.edgeRates[1] = patch.edgeRates[next];
            newPatch1.edgeRates[2] = edgeOppT0;
            newPatch1.edgeRates[3] = patch.edgeRates[prev];

            printf("Split face %d -> newPatch0 corners (uv): (%.4f,%.4f) (%.4f,%.4f) (%.4f,%.4f) "
                   "(%.4f,%.4f) edgeRates: %d %d %d %d\n",
                   patch.ptexFace,
                   newPatch0.parametricCorners[0].x,
                   newPatch0.parametricCorners[0].y,
                   newPatch0.parametricCorners[1].x,
                   newPatch0.parametricCorners[1].y,
                   newPatch0.parametricCorners[2].x,
                   newPatch0.parametricCorners[2].y,
                   newPatch0.parametricCorners[3].x,
                   newPatch0.parametricCorners[3].y,
                   newPatch0.edgeRates[0],
                   newPatch0.edgeRates[1],
                   newPatch0.edgeRates[2],
                   newPatch0.edgeRates[3]);
            printf("Split face %d -> newPatch1 corners (uv): (%.4f,%.4f) (%.4f,%.4f) (%.4f,%.4f) "
                   "(%.4f,%.4f) edgeRates: %d %d %d %d\n",
                   patch.ptexFace,
                   newPatch1.parametricCorners[0].x,
                   newPatch1.parametricCorners[0].y,
                   newPatch1.parametricCorners[1].x,
                   newPatch1.parametricCorners[1].y,
                   newPatch1.parametricCorners[2].x,
                   newPatch1.parametricCorners[2].y,
                   newPatch1.parametricCorners[3].x,
                   newPatch1.parametricCorners[3].y,
                   newPatch1.edgeRates[0],
                   newPatch1.edgeRates[1],
                   newPatch1.edgeRates[2],
                   newPatch1.edgeRates[3]);

            Split(params, newPatch0);
            Split(params, newPatch1);

            break;
        }
    }
}

} // namespace DiagSplit

void Subdivision(Scene *scene, const SubdivisionMesh &mesh, int refineLevel)
{
    using TopologyDescriptor = Far::TopologyDescriptor;

    HostMemoryArena arena;

    const int numFaces = static_cast<int>(mesh.vertsPerFace.size());
    MemoryView<TopologyDescriptor::FVarChannel> channels =
        arena.PushArray<TopologyDescriptor::FVarChannel>(mesh.attributeEnd - mesh.attributeStart);
    size_t numChannels = 0;

    for (size_t i = mesh.attributeStart; i < mesh.attributeEnd; i++)
    {
        Attribute &attribute = scene->attributes[i];
        if (AttributeTypeIsFloat(attribute.type))
        {
            size_t numValues = AttributeTypeGetSize(attribute.type) / sizeof(float);
            channels[numChannels].numValues = numValues;
            channels[numChannels].valueIndices = attribute.indices.data();
            numChannels++;
        }
    }

    TopologyDescriptor desc = {};
    desc.numVertices = static_cast<int>(mesh.vertices.size());
    desc.numFaces = numFaces;
    desc.numVertsPerFace = mesh.vertsPerFace.data();
    desc.vertIndicesPerFace = mesh.indices.data();
#if 0
    desc.numFVarChannels = numChannels;
    desc.fvarChannels = channels.data();
#endif
    desc.numCreases = mesh.creaseSharpnesses.size();
    desc.creaseWeights = mesh.creaseSharpnesses.data();
    // TODO: does the spec match?
    desc.creaseVertexIndexPairs = mesh.creaseIndices.data();
    desc.numCorners = mesh.cornerSharpnesses.size();
    desc.cornerWeights = mesh.cornerSharpnesses.data();
    desc.cornerVertexIndices = mesh.cornerIndices.data();
    desc.numHoles = mesh.holeIndices.size();
    desc.holeIndices = mesh.holeIndices.data();

    Sdc::SchemeType scheme = Sdc::SCHEME_CATMARK;
    Sdc::Options options;
    options.SetVtxBoundaryInterpolation(ToOsdVtxBoundary(mesh.interpolationRule));
    options.SetFVarLinearInterpolation(ToOsdFVarLinear(mesh.fvarLinearInterpolation));
    // options.SetCreasingMethod(CreasingMethod c)
    // options.SetTriangleSubdivision(TriangleSubdivision t)

    Far::TopologyRefiner *refiner = Far::TopologyRefinerFactory<TopologyDescriptor>::Create(
        desc, Far::TopologyRefinerFactory<TopologyDescriptor>::Options(scheme, options));
    YBI_ASSERT(refiner);

    Far::TopologyRefiner::AdaptiveOptions adaptiveOptions(1);
    refiner->RefineAdaptive(adaptiveOptions);

    bool hasCreasesOrCorners = bool(desc.numCreases) || bool(desc.numCorners);
    Far::PatchTableFactory::Options patchOptions(refineLevel);
    patchOptions.endCapType = Far::PatchTableFactory::Options::ENDCAP_GREGORY_BASIS;
    patchOptions.useInfSharpPatch = hasCreasesOrCorners;
    patchOptions.generateFVarLegacyLinearPatches = false;

    const Far::PatchTable *patchTable = Far::PatchTableFactory::Create(*refiner, patchOptions);
    Far::PatchMap patchMap(*patchTable);
    YBI_ASSERT(patchTable);

    Far::PtexIndices ptexIndices(*refiner);
    const int numPtexFaces = ptexIndices.GetNumFaces();
    MemoryView<OsdData<float3>> positions =
        InterpolateVertex<float3>(arena, mesh.vertices, refiner, patchTable);

    if (numPtexFaces == numFaces)
    {
        DiagSplit::DiagSplitParams params = {
            &patchMap,
            patchTable,
            positions,
            3,
            1.f,
            1,
            float4x4::Identity(),
            1024,
            1024,
        };

        for (int face = 0; face < numFaces; face++)
        {
            DiagSplit::SubPatch patch(face);
            DiagSplit::Split(params, patch);
        }
    }

#if 0

    MemoryView<const Attribute *> vertexAttributes =
        arena.PushArray<const Attribute *>(mesh.attributeEnd - mesh.attributeStart);
    size_t numVertexAttributes = 0;
    MemoryView<const Attribute *> varyingAttributes =
        arena.PushArray<const Attribute *>(mesh.attributeEnd - mesh.attributeStart);
    size_t numVaryingAttributes = 0;
    MemoryView<const Attribute *> faceVaryingAttributes =
        arena.PushArray<const Attribute *>(mesh.attributeEnd - mesh.attributeStart);
    size_t numFaceVaryingAttributes = 0;

    auto categorizeAttributes = [&](auto &arr, size_t &count, PrimvarInterpolation interp) {
        for (size_t idx = mesh.attributeStart; idx < mesh.attributeEnd; idx++)
        {
            const Attribute &attr = scene->attributes[idx];
            if (attr.interpolation == interp)
                arr[count++] = &attr;
        }
    };

    categorizeAttributes(vertexAttributes, numVertexAttributes, PrimvarInterpolation::Vertex);
    categorizeAttributes(varyingAttributes, numVaryingAttributes, PrimvarInterpolation::Varying);
    categorizeAttributes(
        faceVaryingAttributes, numFaceVaryingAttributes, PrimvarInterpolation::FaceVarying);

    // CALL_INTERPOLATE(Vertex, InterpolateVertex, vertexAttributes, numVertexAttributes);
    // CALL_INTERPOLATE(Varying, InterpolateVarying, varyingAttributes, numVaryingAttributes);
    // CALL_INTERPOLATE(
    //     FaceVarying, InterpolateFaceVarying, faceVaryingAttributes, numFaceVaryingAttributes);
#endif
}

YBI_NAMESPACE_END
