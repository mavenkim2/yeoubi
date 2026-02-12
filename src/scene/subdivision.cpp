#include "scene/subdivision.h"
#include "scene/attributes.h"
#include "scene/scene.h"
#include "util/base.h"
#include "util/host_memory_arena.h"
#include <memory>
#include <opensubdiv/far/patchTable.h>
#include <opensubdiv/far/patchTableFactory.h>
#include <opensubdiv/far/primvarRefiner.h>
#include <opensubdiv/far/topologyDescriptor.h>
#include <opensubdiv/sdc/types.h>
#include <opensubdiv/version.h>
#include <vector>

YBI_NAMESPACE_BEGIN

namespace Far = OpenSubdiv::Far;
namespace Sdc = OpenSubdiv::Sdc;

template <int numFloats>
struct OsdData
{
    float p[numFloats];

    void Clear()
    {
        for (int i = 0; i < numFloats; i++)
        {
            p[i] = 0.f;
        }
    }
    void AddWithWeight(const OsdData &src, float w)
    {
        for (int i = 0; i < numFloats; i++)
        {
            p[i] += w * src.p[i];
        }
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

template <int numFloats>
static void InterpolateVertex(MemoryView<OsdData<numFloats>> &view,
                              const Far::TopologyRefiner *refiner,
                              const Far::PatchTable *patchTable)
{
    const int numLevels = refiner->GetNumLevels();
    OsdData<numFloats> *src = view.data();

    for (int level = 1; level < numLevels; level++)
    {
        OsdData<numFloats> *dst = src + refiner->GetLevel(level - 1).GetNumVertices();
        Far::PrimvarRefiner(*refiner).Interpolate(level, src, dst);
        src = dst;
    }
    patchTable->ComputeLocalPointValues(src, view.data() + refiner->GetNumVerticesTotal());
}

template <int numFloats>
static void InterpolateVarying(MemoryView<OsdData<numFloats>> &view,
                               const Far::TopologyRefiner *refiner,
                               const Far::PatchTable *patchTable,
                               int channel)
{
    const int numLevels = refiner->GetNumLevels();
    OsdData<numFloats> *src = view.data();

    for (int level = 1; level < numLevels; level++)
    {
        OsdData<numFloats> *dst = src + refiner->GetLevel(level - 1).GetNumVertices();
        Far::PrimvarRefiner(*refiner).InterpolateVarying(level, src, dst);
        src = dst;
    }
    patchTable->ComputeLocalPointValuesVarying(
        src, view.data() + refiner->GetNumFVarValuesTotal(channel));
}

void Subdivision(const SubdivisionMesh &mesh, int refineLevel)
{
    using TopologyDescriptor = Far::TopologyDescriptor;

    HostMemoryArena arena;

    const int numFaces = static_cast<int>(mesh.vertsPerFace.size());

    TopologyDescriptor desc;
    desc.numVertices = static_cast<int>(mesh.vertices.size());
    desc.numFaces = numFaces;
    desc.numFVarChannels = 0;
    desc.fvarChannels = nullptr;

    desc.numVertsPerFace = mesh.vertsPerFace.data();
    desc.vertIndicesPerFace = mesh.indices.data();

    Sdc::SchemeType scheme = Sdc::SCHEME_CATMARK;
    Sdc::Options options;
    options.SetVtxBoundaryInterpolation(ToOsdVtxBoundary(mesh.interpolationRule));
    options.SetFVarLinearInterpolation(ToOsdFVarLinear(mesh.fvarLinearInterpolation));

    Far::TopologyRefiner *refiner = Far::TopologyRefinerFactory<TopologyDescriptor>::Create(
        desc, Far::TopologyRefinerFactory<TopologyDescriptor>::Options(scheme, options));
    YBI_ASSERT(refiner);

    Far::TopologyRefiner::AdaptiveOptions adaptiveOptions(refineLevel);
    refiner->RefineAdaptive(adaptiveOptions);

    bool hasCreasesOrCorners = mesh.cornerIndices.size() || mesh.creaseIndices.size();
    Far::PatchTableFactory::Options patchOptions(refineLevel);
    patchOptions.endCapType = Far::PatchTableFactory::Options::ENDCAP_GREGORY_BASIS;
    patchOptions.useInfSharpPatch = true;
    patchOptions.generateFVarLegacyLinearPatches = false;

    const Far::PatchTable *patchTable = Far::PatchTableFactory::Create(*refiner, patchOptions);
    YBI_ASSERT(patchTable);

    const int numVerticesTotal = refiner->GetNumVerticesTotal();
    const int numLocalPointsTotal = patchTable->GetNumLocalPoints();

    MemoryView<OsdData<3>> positions =
        arena.PushArray<OsdData<3>>(numVerticesTotal + numLocalPointsTotal);
    const int numCoarse = refiner->GetLevel(0).GetNumVertices();
    YBI_ASSERT(numCoarse == mesh.vertices.size());

    memcpy(positions.data(), mesh.vertices.data(), sizeof(float3) * mesh.vertices.size());
    InterpolateVertex<3>(positions, refiner, patchTable);

#if 0


    pt->ComputeLocalPointValues(src, verts.data() + numVerticesTotal);

    const size_t numPatchesTotal = pt->GetNumPatchesTotal();
    result.limitPositions.Resize(numPatchesTotal);

    const int numArrays = pt->GetNumPatchArrays();
    int globalPatchIndex = 0;
    std::vector<float> wP;

    for (int array = 0; array < numArrays; array++)
    {
        const int numPatchesInArray = pt->GetNumPatches(array);
        const int numCVsPerPatch = pt->GetPatchArrayDescriptor(array).GetNumControlVertices();
        wP.resize(numCVsPerPatch);

        for (int patch = 0; patch < numPatchesInArray; patch++)
        {
            Far::PatchTable::PatchHandle handle;
            handle.arrayIndex = array;
            handle.patchIndex = globalPatchIndex;
            handle.vertIndex = patch * numCVsPerPatch;

            pt->EvaluateBasis(handle, 0.5f, 0.5f, wP.data());

            Far::ConstIndexArray cvs = pt->GetPatchVertices(array, patch);
            float x = 0.f, y = 0.f, z = 0.f;
            for (int cv = 0; cv < cvs.size(); cv++)
            {
                int idx = cvs[cv];
                float w = wP[cv];
                const OsdVertex &v = verts[idx];
                x += w * v.p[0];
                y += w * v.p[1];
                z += w * v.p[2];
            }
            result.limitPositions[globalPatchIndex] = make_float3(x, y, z);
            ++globalPatchIndex;
        }
    }

#endif
}

YBI_NAMESPACE_END
