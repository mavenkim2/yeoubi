#include "io/usd_subdiv_json_io.h"

#include <opensubdiv/far/patchTableFactory.h>
#include <opensubdiv/far/topologyDescriptor.h>
#include <opensubdiv/far/topologyRefinerFactory.h>

#include <algorithm>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <string>
#include <unordered_map>
#include <vector>

using namespace OpenSubdiv;

struct TriMesh
{
    pxr::VtVec3fArray positions;
    std::vector<int> indices;
};

struct CreasePairs
{
    std::vector<int> pairs;
    std::vector<float> weights;
};

static CreasePairs BuildCreasePairs(const SelectedSubdivMesh &m)
{
    CreasePairs out = {};
    size_t c = 0;
    for (size_t i = 0; i < m.creaseLengths.size(); i++)
    {
        const int len = m.creaseLengths[i];
        if (len < 2 || c + len > m.creaseIndices.size())
        {
            c += std::max(0, len);
            continue;
        }
        const float w = (i < m.creaseSharpnesses.size()) ? m.creaseSharpnesses[i] : 0.0f;
        for (int j = 0; j + 1 < len; j++)
        {
            out.pairs.push_back(m.creaseIndices[c + j]);
            out.pairs.push_back(m.creaseIndices[c + j + 1]);
            out.weights.push_back(w);
        }
        c += len;
    }
    return out;
}

static Sdc::SchemeType SchemeFromString(const std::string &s)
{
    if (s == "loop")
    {
        return Sdc::SCHEME_LOOP;
    }
    if (s == "bilinear")
    {
        return Sdc::SCHEME_BILINEAR;
    }
    return Sdc::SCHEME_CATMARK;
}

static Sdc::Options::VtxBoundaryInterpolation VtxBoundaryFromString(const std::string &s)
{
    if (s == "none")
    {
        return Sdc::Options::VTX_BOUNDARY_NONE;
    }
    if (s == "edgeOnly")
    {
        return Sdc::Options::VTX_BOUNDARY_EDGE_ONLY;
    }
    return Sdc::Options::VTX_BOUNDARY_EDGE_AND_CORNER;
}

static Sdc::Options::FVarLinearInterpolation FVarLinearFromString(const std::string &s)
{
    if (s == "none")
    {
        return Sdc::Options::FVAR_LINEAR_NONE;
    }
    if (s == "cornersOnly")
    {
        return Sdc::Options::FVAR_LINEAR_CORNERS_ONLY;
    }
    if (s == "cornersPlus1")
    {
        return Sdc::Options::FVAR_LINEAR_CORNERS_PLUS1;
    }
    if (s == "cornersPlus2")
    {
        return Sdc::Options::FVAR_LINEAR_CORNERS_PLUS2;
    }
    if (s == "boundaries")
    {
        return Sdc::Options::FVAR_LINEAR_BOUNDARIES;
    }
    if (s == "all" || s == "bilinear")
    {
        return Sdc::Options::FVAR_LINEAR_ALL;
    }
    return Sdc::Options::FVAR_LINEAR_CORNERS_PLUS1;
}

static Sdc::Options::CreasingMethod CreasingMethodFromString(const std::string &s)
{
    if (s == "chaikin")
    {
        return Sdc::Options::CREASE_CHAIKIN;
    }
    return Sdc::Options::CREASE_UNIFORM;
}

static Sdc::Options::TriangleSubdivision TriangleSubFromString(const std::string &s)
{
    if (s == "smooth")
    {
        return Sdc::Options::TRI_SUB_SMOOTH;
    }
    return Sdc::Options::TRI_SUB_CATMARK;
}

static bool WriteTriObj(const std::string &path, const TriMesh &mesh)
{
    std::ofstream out(path);
    if (!out.is_open())
    {
        return false;
    }
    for (const pxr::GfVec3f &p : mesh.positions)
    {
        out << "v " << p[0] << " " << p[1] << " " << p[2] << "\n";
    }
    for (int i = 0; i + 2 < int(mesh.indices.size()); i += 3)
    {
        out << "f " << (mesh.indices[i + 0] + 1) << " " << (mesh.indices[i + 1] + 1) << " "
            << (mesh.indices[i + 2] + 1) << "\n";
    }
    return true;
}

static TriMesh BuildTriangulatedControlCage(const SelectedSubdivMesh &m)
{
    TriMesh mesh = {};
    mesh.positions = m.points;

    size_t cursor = 0;
    for (size_t f = 0; f < m.faceVertexCounts.size(); ++f)
    {
        const int n = m.faceVertexCounts[f];
        if (n < 3 || cursor + n > m.faceVertexIndices.size())
        {
            cursor += std::max(0, n);
            continue;
        }

        const int i0 = m.faceVertexIndices[cursor + 0];
        for (int i = 1; i + 1 < n; ++i)
        {
            mesh.indices.push_back(i0);
            mesh.indices.push_back(m.faceVertexIndices[cursor + i]);
            mesh.indices.push_back(m.faceVertexIndices[cursor + i + 1]);
        }
        cursor += n;
    }

    return mesh;
}

static int CountControlCageEdgesWithOver2FacesFromIndices(const SelectedSubdivMesh &m)
{
    std::unordered_map<uint64_t, int> edgeUseCount;
    edgeUseCount.reserve(m.faceVertexIndices.size());

    size_t cursor = 0;
    for (size_t f = 0; f < m.faceVertexCounts.size(); ++f)
    {
        const int n = m.faceVertexCounts[f];
        if (n < 2 || cursor + n > m.faceVertexIndices.size())
        {
            cursor += std::max(0, n);
            continue;
        }
        for (int i = 0; i < n; ++i)
        {
            const int a = m.faceVertexIndices[cursor + i];
            const int b = m.faceVertexIndices[cursor + ((i + 1) % n)];
            const int lo = std::min(a, b);
            const int hi = std::max(a, b);
            const uint64_t key = (uint64_t(uint32_t(lo)) << 32) | uint64_t(uint32_t(hi));
            edgeUseCount[key]++;
        }
        cursor += n;
    }

    int over2 = 0;
    for (const auto &it : edgeUseCount)
    {
        if (it.second > 2)
        {
            over2++;
        }
    }
    return over2;
}

int main(int argc, char **argv)
{
    if (argc < 2)
    {
        std::fprintf(stderr, "Usage: %s <selected-subdiv.json> [level>=1] [out.obj]\n", argv[0]);
        return 2;
    }

    const std::string inJson = argv[1];
    int level = 1;
    std::string outObj;
    if (argc >= 3)
    {
        level = std::max(1, std::atoi(argv[2]));
    }
    if (argc >= 4)
    {
        outObj = argv[3];
    }
    if (outObj.empty())
    {
        outObj = "tests/bvh/out/refined_adaptive_level" + std::to_string(level) + ".obj";
    }

    SelectedSubdivMesh m = {};
    UsdCameraInfo camera = {};
    if (!ybi::testio::LoadSelectedSubdivFromJson(inJson, m, camera))
    {
        std::fprintf(stderr, "Failed to load JSON: %s\n", inJson.c_str());
        return 1;
    }

    CreasePairs creases = BuildCreasePairs(m);
    Far::TopologyDescriptor d = {};
    d.numVertices = int(m.points.size());
    d.numFaces = int(m.faceVertexCounts.size());
    d.numVertsPerFace = m.faceVertexCounts.data();
    d.vertIndicesPerFace = m.faceVertexIndices.data();
    d.numCreases = int(creases.weights.size());
    d.creaseVertexIndexPairs = creases.pairs.data();
    d.creaseWeights = creases.weights.data();
    d.numCorners = int(m.cornerIndices.size());
    d.cornerVertexIndices = m.cornerIndices.data();
    d.cornerWeights = m.cornerSharpnesses.data();
    d.numHoles = int(m.holeIndices.size());
    d.holeIndices = m.holeIndices.data();

    Sdc::Options sdcOptions;
    sdcOptions.SetVtxBoundaryInterpolation(VtxBoundaryFromString(m.vertexBoundaryInterpolation));
    sdcOptions.SetFVarLinearInterpolation(FVarLinearFromString(m.fvarLinearInterpolation));
    sdcOptions.SetCreasingMethod(CreasingMethodFromString(m.creasingMethod));
    sdcOptions.SetTriangleSubdivision(TriangleSubFromString(m.triangleSubdivision));

    Far::TopologyRefinerFactory<Far::TopologyDescriptor>::Options o(
        SchemeFromString(m.subdivisionScheme), sdcOptions);

    const int edgesWithOver2Faces = CountControlCageEdgesWithOver2FacesFromIndices(m);
    if (edgesWithOver2Faces > 0)
    {
        std::fprintf(
            stderr,
            "WARNING: non-manifold control cage detected from indices (edges with >2 faces: %d)\n",
            edgesWithOver2Faces);
    }

    Far::TopologyRefiner *refiner =
        Far::TopologyRefinerFactory<Far::TopologyDescriptor>::Create(d, o);
    if (!refiner)
    {
        std::fprintf(stderr, "Failed to create TopologyRefiner\n");
        return 1;
    }

    Far::TopologyRefiner::AdaptiveOptions adaptiveOptions(level);
    refiner->RefineAdaptive(adaptiveOptions);
    if (refiner->GetMaxLevel() < level)
    {
        std::fprintf(stderr, "Adaptive refine did not produce level %d\n", level);
        delete refiner;
        return 1;
    }

    Far::PatchTableFactory::Options patchOptions(level);
    patchOptions.endCapType = Far::PatchTableFactory::Options::ENDCAP_GREGORY_BASIS;
    patchOptions.useInfSharpPatch = (d.numCreases > 0) || (d.numCorners > 0);
    patchOptions.generateFVarLegacyLinearPatches = false;
    const Far::PatchTable *patchTable = Far::PatchTableFactory::Create(*refiner, patchOptions);
    if (!patchTable)
    {
        std::fprintf(stderr, "Failed to create PatchTable\n");
        delete refiner;
        return 1;
    }

    // TODO: restart adaptive tessellation path from patch params; current fallback writes
    // triangulated control cage while keeping refiner/patch table creation validated.
    const TriMesh triMesh = BuildTriangulatedControlCage(m);
    if (!WriteTriObj(outObj, triMesh))
    {
        std::fprintf(stderr, "Failed to write OBJ: %s\n", outObj.c_str());
        delete patchTable;
        delete refiner;
        return 1;
    }

    std::printf("Wrote fallback adaptive OBJ (control cage triangulation): %s\n", outObj.c_str());
    std::printf("  levelRequested=%d\n", level);
    std::printf("  fallback verts=%zu tris=%zu\n",
                triMesh.positions.size(),
                triMesh.indices.size() / 3);
    std::printf("  refinedMaxLevel=%d patches=%d\n",
                refiner->GetMaxLevel(),
                patchTable->GetNumPatchesTotal());
    std::printf("  controlCageEdgesWithOver2Faces=%d\n", edgesWithOver2Faces);

    delete patchTable;
    delete refiner;
    return 0;
}
