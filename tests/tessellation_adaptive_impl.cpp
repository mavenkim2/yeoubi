#include "io/usd_subdiv_json_io.h"

#include <opensubdiv/far/primvarRefiner.h>
#include <opensubdiv/far/ptexIndices.h>
#include <opensubdiv/far/topologyDescriptor.h>
#include <opensubdiv/far/topologyRefinerFactory.h>

#include "util/assert.h"
#include <algorithm>
#include <array>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <string>
#include <unordered_map>
#include <vector>

using namespace OpenSubdiv;

struct CreasePairs
{
    std::vector<int> pairs;
    std::vector<float> weights;
};

struct Primvar3
{
    float x = 0.0f;
    float y = 0.0f;
    float z = 0.0f;
    Primvar3() = default;
    explicit Primvar3(const pxr::GfVec3f &p) : x(p[0]), y(p[1]), z(p[2]) {}
    void Clear()
    {
        x = y = z = 0.0f;
    }
    void AddWithWeight(const Primvar3 &p, float w)
    {
        x += p.x * w;
        y += p.y * w;
        z += p.z * w;
    }
    pxr::GfVec3f ToGf() const
    {
        return pxr::GfVec3f(x, y, z);
    }
};

struct PtexFaceAdj
{
    std::array<int, 4> adjFace = {-1, -1, -1, -1};
    std::array<int, 4> adjEdge = {-1, -1, -1, -1};
    std::array<int, 4> edgeIndex = {-1, -1, -1, -1};
    bool fromNgon = false;
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

static bool WriteObjLevel(const std::string &path,
                          const Far::TopologyLevel &level,
                          const std::vector<Primvar3> &pos,
                          int &nonQuadFacesOut)
{
    std::ofstream out(path);
    if (!out.is_open())
    {
        return false;
    }
    for (const auto &p : pos)
    {
        const pxr::GfVec3f g = p.ToGf();
        out << "v " << g[0] << " " << g[1] << " " << g[2] << "\n";
    }
    nonQuadFacesOut = 0;
    const int numVerts = int(pos.size());
    for (int f = 0; f < level.GetNumFaces(); f++)
    {
        auto verts = level.GetFaceVertices(f);
        if (verts.size() != 4)
        {
            nonQuadFacesOut++;
        }
        for (int i = 0; i < verts.size(); i++)
        {
            if (verts[i] < 0 || verts[i] >= numVerts)
            {
                std::fprintf(stderr,
                             "Invalid face index: face=%d corner=%d index=%d numVerts=%d\n",
                             f,
                             i,
                             verts[i],
                             numVerts);
                return false;
            }
        }
        out << "f";
        for (int i = 0; i < verts.size(); i++)
        {
            out << " " << (verts[i] + 1);
        }
        out << "\n";
    }
    return true;
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

static int BuildUniquePtexEdgeIds(const Far::TopologyRefiner &refiner,
                                  const Far::TopologyLevel &level0,
                                  const Far::PtexIndices &ptex,
                                  int &numPtexFacesOut)
{
    numPtexFacesOut = ptex.GetNumFaces();
    std::vector<PtexFaceAdj> faces(std::max(0, numPtexFacesOut));

    for (int f = 0; f < level0.GetNumFaces(); ++f)
    {
        const int n = level0.GetFaceVertices(f).size();
        const int qCount = (n == 4) ? 1 : n;
        const int basePtex = ptex.GetFaceId(f);
        for (int q = 0; q < qCount; ++q)
        {
            const int pf = basePtex + q;
            if (pf < 0 || pf >= numPtexFacesOut)
            {
                continue;
            }
            int adjFaces[4] = {-1, -1, -1, -1};
            int adjEdges[4] = {-1, -1, -1, -1};
            ptex.GetAdjacency(refiner, f, q, adjFaces, adjEdges);
            for (int e = 0; e < 4; ++e)
            {
                faces[pf].adjFace[e] = adjFaces[e];
                faces[pf].adjEdge[e] = adjEdges[e];
            }
            faces[pf].fromNgon = (n != 4);
        }
    }

    int edgeCount = 0;
    for (int pf = 0; pf < numPtexFacesOut; ++pf)
    {
        for (int e = 0; e < 4; ++e)
        {
            const int af = faces[pf].adjFace[e];
            const int ae = faces[pf].adjEdge[e];
            YBI_ASSERT(af < numPtexFacesOut && ae >= 0 && ae < 4);

            if (faces[pf].edgeIndex[e] >= 0)
            {
                continue;
            }
            const int neighborEdge = faces[af].edgeIndex[ae];
            if (neighborEdge >= 0)
            {
                faces[pf].edgeIndex[e] = neighborEdge;
                continue;
            }
            const int edgeId = edgeCount++;
            faces[pf].edgeIndex[e] = edgeId;
            faces[af].edgeIndex[ae] = edgeId;
        }
    }
    return edgeCount;
}

int main(int argc, char **argv)
{
    if (argc < 2 || argc > 4)
    {
        std::fprintf(stderr, "Usage: %s <selected-subdiv.json> [level>=1] [out.obj]\n", argv[0]);
        return 2;
    }

    const std::string inJson = argv[1];
    const int level = (argc >= 3) ? std::max(1, std::atoi(argv[2])) : 1;
    const std::string outObj =
        (argc >= 4) ? argv[3]
                    : ("tests/bvh/out/refined_uniform_level" + std::to_string(level) + ".obj");

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

    Far::TopologyRefinerFactory<Far::TopologyDescriptor>::Options o(
        SchemeFromString(m.subdivisionScheme));

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

    const Far::TopologyLevel &level0 = refiner->GetLevel(0);

    refiner->RefineUniform(Far::TopologyRefiner::UniformOptions(level));
    if (refiner->GetMaxLevel() < level)
    {
        std::fprintf(stderr, "Uniform refine did not produce level %d\n", level);
        delete refiner;
        return 1;
    }

    std::vector<Primvar3> p0;
    p0.reserve(m.points.size());
    for (const auto &p : m.points)
    {
        p0.emplace_back(p);
    }
    std::vector<Primvar3> p1(refiner->GetLevel(level).GetNumVertices());
    Far::PrimvarRefiner prim(*refiner);
    prim.Interpolate(level, p0, p1);

    Far::PtexIndices ptex(*refiner);
    int ptexFaceCount = 0;
    BuildUniquePtexEdgeIds(*refiner, level0, ptex, ptexFaceCount);

    int nonQuadFaces = 0;
    if (!WriteObjLevel(outObj, refiner->GetLevel(level), p1, nonQuadFaces))
    {
        std::fprintf(stderr, "Failed to write OBJ: %s\n", outObj.c_str());
        delete refiner;
        return 1;
    }

    std::printf("Wrote uniform level-%d OBJ: %s\n", level, outObj.c_str());
    std::printf("  verts=%zu faces=%d nonQuads=%d\n",
                p1.size(),
                refiner->GetLevel(level).GetNumFaces(),
                nonQuadFaces);
    std::printf("  ptexFaces=%d\n", ptexFaceCount);
    std::printf("  controlCageEdgesWithOver2Faces=%d\n", edgesWithOver2Faces);

    delete refiner;
    return 0;
}
