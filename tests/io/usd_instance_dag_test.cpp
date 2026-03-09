#include "io/usd/instance_dag_build.h"
#include "io/usd/load.h"
#include "scene/scene.h"

#include <pxr/base/gf/matrix4d.h>
#include <pxr/base/gf/vec3d.h>
#include <pxr/usd/sdf/layer.h>
#include <pxr/usd/usd/stage.h>

#include <cmath>
#include <cstdio>
#include <fstream>
#include <string>
#include <unordered_map>

namespace
{

static const char *kNoInstancesFixture = R"USDA(#usda 1.0
def Xform "GeomRoot"
{
    def Mesh "MeshA"
    {
        int[] faceVertexCounts = [3]
        int[] faceVertexIndices = [0, 1, 2]
        point3f[] points = [(0, 0, 0), (1, 0, 0), (0, 1, 0)]
    }
}
)USDA";

static const char *kNormalInstancingFixture = R"USDA(#usda 1.0
def Xform "Model"
{
    def Mesh "MeshA"
    {
        int[] faceVertexCounts = [3]
        int[] faceVertexIndices = [0, 1, 2]
        point3f[] points = [(0, 0, 0), (1, 0, 0), (0, 1, 0)]
    }
}

def Xform "InstA" (
    instanceable = true
    prepend references = </Model>
)
{
}

def Xform "InstB" (
    instanceable = true
    prepend references = </Model>
)
{
}
)USDA";

static const char *kNestedInstancingFixture = R"USDA(#usda 1.0
def Xform "LeafModel"
{
    def Mesh "LeafMesh"
    {
        int[] faceVertexCounts = [3]
        int[] faceVertexIndices = [0, 1, 2]
        point3f[] points = [(0, 0, 0), (1, 0, 0), (0, 1, 0)]
    }
}

def Xform "PIProtoWithInstance"
{
    def Xform "NestedInst" (
        instanceable = true
        prepend references = </LeafModel>
    )
    {
    }
}

def Xform "PIProtoPlain"
{
    def Mesh "PlainMesh"
    {
        int[] faceVertexCounts = [3]
        int[] faceVertexIndices = [0, 1, 2]
        point3f[] points = [(0, 0, 0), (0, 1, 0), (0, 0, 1)]
    }
}

def Xform "TopModel"
{
    def Mesh "TopMesh"
    {
        int[] faceVertexCounts = [3]
        int[] faceVertexIndices = [0, 1, 2]
        point3f[] points = [(0, 0, 0), (1, 0, 0), (0, 0, 1)]
    }

    def PointInstancer "PI"
    {
        rel prototypes = [</PIProtoWithInstance>, </PIProtoPlain>]
        int[] protoIndices = [0, 1]
        point3f[] positions = [(0, 0, 0), (3, 0, 0)]
    }
}

def Xform "TopInst" (
    instanceable = true
    prepend references = </TopModel>
)
{
}
)USDA";

static const char *kPointInstancerTransformFixture = R"USDA(#usda 1.0
def Xform "World"
{
    double3 xformOp:translate = (10, 0, 0)
    uniform token[] xformOpOrder = ["xformOp:translate"]

    def PointInstancer "PI"
    {
        def Scope "Prototypes"
        {
            def Xform "ProtoA"
            {
                double3 xformOp:translate = (0, 5, 0)
                uniform token[] xformOpOrder = ["xformOp:translate"]
            }
        }

        rel prototypes = [</World/PI/Prototypes/ProtoA>]
        int[] protoIndices = [0, 0]
        point3f[] positions = [(1, 0, 0), (0, 2, 0)]
        float3[] scales = [(2, 2, 2), (1, 3, 1)]
    }
}
)USDA";

static const char *kVisibilityFixture = R"USDA(#usda 1.0
def Xform "World"
{
    def Mesh "VisibleMesh"
    {
        int[] faceVertexCounts = [3]
        int[] faceVertexIndices = [0, 1, 2]
        point3f[] points = [(0, 0, 0), (1, 0, 0), (0, 1, 0)]
        uniform token subdivisionScheme = "none"
    }

    def Mesh "HiddenMesh"
    {
        token visibility = "invisible"
        int[] faceVertexCounts = [3]
        int[] faceVertexIndices = [0, 1, 2]
        point3f[] points = [(0, 0, 1), (1, 0, 1), (0, 1, 1)]
        uniform token subdivisionScheme = "none"
    }

    def Xform "HiddenGroup"
    {
        token visibility = "invisible"

        def Mesh "HiddenChildMesh"
        {
            int[] faceVertexCounts = [3]
            int[] faceVertexIndices = [0, 1, 2]
            point3f[] points = [(0, 0, 2), (1, 0, 2), (0, 1, 2)]
            uniform token subdivisionScheme = "none"
        }
    }
}
)USDA";

static pxr::UsdStageRefPtr CreateStageFromUsda(const char *usda, std::string *error)
{
    pxr::UsdStageRefPtr stage = pxr::UsdStage::CreateInMemory();
    if (!stage)
    {
        if (error)
        {
            *error = "failed to create in-memory stage";
        }
        return {};
    }

    const pxr::SdfLayerRefPtr layer = stage->GetRootLayer();
    if (!layer || !layer->ImportFromString(usda))
    {
        if (error)
        {
            *error = "failed to import fixture USDA";
        }
        return {};
    }

    return stage;
}

static ybi::Float3x4 ConvertAffineTransformForTest(const pxr::GfMatrix4d &m)
{
    const pxr::GfMatrix4d t = m.GetTranspose();
    const pxr::GfVec4d r0 = t.GetRow(0);
    const pxr::GfVec4d r1 = t.GetRow(1);
    const pxr::GfVec4d r2 = t.GetRow(2);
    return ybi::Float3x4(static_cast<float>(r0[0]),
                         static_cast<float>(r0[1]),
                         static_cast<float>(r0[2]),
                         static_cast<float>(r0[3]),
                         static_cast<float>(r1[0]),
                         static_cast<float>(r1[1]),
                         static_cast<float>(r1[2]),
                         static_cast<float>(r1[3]),
                         static_cast<float>(r2[0]),
                         static_cast<float>(r2[1]),
                         static_cast<float>(r2[2]),
                         static_cast<float>(r2[3]));
}

static bool MatrixNear(const ybi::Float3x4 &actual,
                       const ybi::Float3x4 &expected,
                       float eps,
                       const char *label)
{
    for (int r = 0; r < 3; r++)
    {
        for (int c = 0; c < 4; c++)
        {
            if (std::fabs(actual.m[r][c] - expected.m[r][c]) > eps)
            {
                std::fprintf(stderr,
                             "%s: matrix mismatch at [%d][%d], got %.6f expected %.6f\n",
                             label,
                             r,
                             c,
                             actual.m[r][c],
                             expected.m[r][c]);
                return false;
            }
        }
    }
    return true;
}

static bool CheckCommonDAGInvariants(const ybi::USDBuildSceneDAG &dag, std::string *error)
{
    if (dag.rootSceneIndex >= dag.scenes.size())
    {
        if (error)
        {
            *error = "rootSceneIndex out of range";
        }
        return false;
    }

    for (size_t sceneIndex = 0; sceneIndex < dag.scenes.size(); sceneIndex++)
    {
        const ybi::USDBuildScene &scene = dag.scenes[sceneIndex];
        for (size_t i = 0; i < scene.instances.size(); i++)
        {
            if (scene.instances[i].childSceneIndex >= dag.scenes.size())
            {
                if (error)
                {
                    *error = "instance childSceneIndex out of range in scene " +
                             std::to_string(sceneIndex);
                }
                return false;
            }
        }
    }
    return true;
}

static bool RunPointInstancerTransformLoadTest()
{
    const std::string fixturePath = "tests/io/_tmp_point_instancer_transform_fixture.usda";
    {
        std::ofstream fixtureFile(fixturePath, std::ios::out | std::ios::binary);
        if (!fixtureFile.is_open())
        {
            std::fprintf(stderr, "PointInstancerTransforms: failed to open fixture file\n");
            return false;
        }
        fixtureFile << kPointInstancerTransformFixture;
    }

    ybi::ScenePool scenePool = {};
    ybi::LoadUSDScene(&scenePool, fixturePath);
    std::remove(fixturePath.c_str());

    if (scenePool.scenes.empty() || scenePool.rootSceneIndex >= scenePool.scenes.size())
    {
        std::fprintf(stderr, "PointInstancerTransforms: invalid scene pool after LoadUSDScene\n");
        return false;
    }

    const ybi::Scene &root = *scenePool.scenes[scenePool.rootSceneIndex];
    if (root.instances.size() != 2)
    {
        std::fprintf(stderr,
                     "PointInstancerTransforms: expected 2 root instances, got %zu\n",
                     root.instances.size());
        return false;
    }
    if (root.instances[0].childSceneIndex >= root.childScenes.size() ||
        root.instances[1].childSceneIndex >= root.childScenes.size())
    {
        std::fprintf(stderr,
                     "PointInstancerTransforms: childSceneIndex out of range in root scene\n");
        return false;
    }
    if (root.childScenes[root.instances[0].childSceneIndex] !=
        root.childScenes[root.instances[1].childSceneIndex])
    {
        std::fprintf(stderr,
                     "PointInstancerTransforms: expected both instances to reference same child "
                     "scene pointer\n");
        return false;
    }

    pxr::GfMatrix4d pointInstancerLocalToWorld(1.0);
    pointInstancerLocalToWorld.SetTranslateOnly(pxr::GfVec3d(10.0, 0.0, 0.0));

    pxr::GfMatrix4d prototypeLocalToParent(1.0);
    prototypeLocalToParent.SetTranslateOnly(pxr::GfVec3d(0.0, 5.0, 0.0));

    pxr::GfMatrix4d t0(1.0);
    t0.SetTranslateOnly(pxr::GfVec3d(1.0, 0.0, 0.0));
    pxr::GfMatrix4d s0(1.0);
    s0.SetScale(pxr::GfVec3d(2.0, 2.0, 2.0));
    const ybi::Float3x4 expected0 = ConvertAffineTransformForTest(
        pointInstancerLocalToWorld * t0 * s0 * prototypeLocalToParent);

    pxr::GfMatrix4d t1(1.0);
    t1.SetTranslateOnly(pxr::GfVec3d(0.0, 2.0, 0.0));
    pxr::GfMatrix4d s1(1.0);
    s1.SetScale(pxr::GfVec3d(1.0, 3.0, 1.0));
    const ybi::Float3x4 expected1 = ConvertAffineTransformForTest(
        pointInstancerLocalToWorld * t1 * s1 * prototypeLocalToParent);

    if (!MatrixNear(root.instances[0].parentFromLocal, expected0, 1e-4f, "PointInstancer[0]"))
    {
        return false;
    }
    if (!MatrixNear(root.instances[1].parentFromLocal, expected1, 1e-4f, "PointInstancer[1]"))
    {
        return false;
    }

    return true;
}

static bool RunNoInstancesTest()
{
    std::string error;
    pxr::UsdStageRefPtr stage = CreateStageFromUsda(kNoInstancesFixture, &error);
    if (!stage)
    {
        std::fprintf(stderr, "NoInstances: %s\n", error.c_str());
        return false;
    }

    ybi::USDBuildSceneDAG dag = {};
    if (!ybi::BuildInstanceDAGFromUSD(stage, &dag, &error))
    {
        std::fprintf(stderr, "NoInstances: BuildInstanceDAGFromUSD failed: %s\n", error.c_str());
        return false;
    }
    if (!CheckCommonDAGInvariants(dag, &error))
    {
        std::fprintf(stderr, "NoInstances: invariant failed: %s\n", error.c_str());
        return false;
    }

    if (dag.scenes.size() != 1)
    {
        std::fprintf(stderr, "NoInstances: expected 1 scene, got %zu\n", dag.scenes.size());
        return false;
    }

    const ybi::USDBuildScene &root = dag.scenes[dag.rootSceneIndex];
    if (!root.instances.empty())
    {
        std::fprintf(stderr, "NoInstances: expected 0 instances in root, got %zu\n",
                     root.instances.size());
        return false;
    }
    if (root.meshes.empty())
    {
        std::fprintf(stderr, "NoInstances: expected at least one mesh in root\n");
        return false;
    }

    return true;
}

static bool RunVisibilityCullTest()
{
    std::string error;
    pxr::UsdStageRefPtr stage = CreateStageFromUsda(kVisibilityFixture, &error);
    if (!stage)
    {
        std::fprintf(stderr, "VisibilityCull: %s\n", error.c_str());
        return false;
    }

    ybi::USDBuildSceneDAG dag = {};
    if (!ybi::BuildInstanceDAGFromUSD(stage, &dag, &error))
    {
        std::fprintf(stderr, "VisibilityCull: BuildInstanceDAGFromUSD failed: %s\n", error.c_str());
        return false;
    }
    if (!CheckCommonDAGInvariants(dag, &error))
    {
        std::fprintf(stderr, "VisibilityCull: invariant failed: %s\n", error.c_str());
        return false;
    }

    const ybi::USDBuildScene &root = dag.scenes[dag.rootSceneIndex];
    if (root.meshes.size() != 1)
    {
        std::fprintf(stderr, "VisibilityCull: expected 1 visible mesh in DAG, got %zu\n",
                     root.meshes.size());
        return false;
    }
    if (root.meshes[0].path.find("VisibleMesh") == std::string::npos)
    {
        std::fprintf(stderr, "VisibilityCull: wrong visible mesh path: %s\n",
                     root.meshes[0].path.c_str());
        return false;
    }

    const std::string fixturePath = "tests/io/_tmp_visibility_fixture.usda";
    {
        std::ofstream fixtureFile(fixturePath, std::ios::out | std::ios::binary);
        if (!fixtureFile.is_open())
        {
            std::fprintf(stderr, "VisibilityCull: failed to open fixture file\n");
            return false;
        }
        fixtureFile << kVisibilityFixture;
    }

    ybi::ScenePool scenePool = {};
    ybi::LoadUSDScene(&scenePool, fixturePath);
    std::remove(fixturePath.c_str());

    if (scenePool.scenes.empty() || scenePool.rootSceneIndex >= scenePool.scenes.size())
    {
        std::fprintf(stderr, "VisibilityCull: invalid scene pool after LoadUSDScene\n");
        return false;
    }

    const ybi::Scene &runtimeRoot = *scenePool.scenes[scenePool.rootSceneIndex];
    if (runtimeRoot.meshes.size() != 1)
    {
        std::fprintf(stderr, "VisibilityCull: expected 1 visible mesh at runtime, got %zu\n",
                     runtimeRoot.meshes.size());
        return false;
    }

    return true;
}

static bool RunNormalInstancingTest()
{
    std::string error;
    pxr::UsdStageRefPtr stage = CreateStageFromUsda(kNormalInstancingFixture, &error);
    if (!stage)
    {
        std::fprintf(stderr, "NormalInstancing: %s\n", error.c_str());
        return false;
    }

    ybi::USDBuildSceneDAG dag = {};
    if (!ybi::BuildInstanceDAGFromUSD(stage, &dag, &error))
    {
        std::fprintf(stderr, "NormalInstancing: BuildInstanceDAGFromUSD failed: %s\n",
                     error.c_str());
        return false;
    }
    if (!CheckCommonDAGInvariants(dag, &error))
    {
        std::fprintf(stderr, "NormalInstancing: invariant failed: %s\n", error.c_str());
        return false;
    }

    if (dag.scenes.size() < 2)
    {
        std::fprintf(stderr, "NormalInstancing: expected at least 2 scenes, got %zu\n",
                     dag.scenes.size());
        return false;
    }

    const ybi::USDBuildScene &root = dag.scenes[dag.rootSceneIndex];
    if (root.instances.size() != 2)
    {
        std::fprintf(stderr, "NormalInstancing: expected 2 root instances, got %zu\n",
                     root.instances.size());
        return false;
    }

    if (root.instances[0].childSceneIndex != root.instances[1].childSceneIndex)
    {
        std::fprintf(stderr, "NormalInstancing: expected both instances to share child scene\n");
        return false;
    }

    return true;
}

static bool RunNestedInstancingTest()
{
    std::string error;
    pxr::UsdStageRefPtr stage = CreateStageFromUsda(kNestedInstancingFixture, &error);
    if (!stage)
    {
        std::fprintf(stderr, "NestedInstancing: %s\n", error.c_str());
        return false;
    }

    ybi::USDBuildSceneDAG dag = {};
    if (!ybi::BuildInstanceDAGFromUSD(stage, &dag, &error))
    {
        std::fprintf(stderr, "NestedInstancing: BuildInstanceDAGFromUSD failed: %s\n",
                     error.c_str());
        return false;
    }
    if (!CheckCommonDAGInvariants(dag, &error))
    {
        std::fprintf(stderr, "NestedInstancing: invariant failed: %s\n", error.c_str());
        return false;
    }

    const ybi::USDBuildScene &root = dag.scenes[dag.rootSceneIndex];
    if (root.instances.empty())
    {
        std::fprintf(stderr, "NestedInstancing: expected root to contain instances\n");
        return false;
    }

    bool foundPrototypeWithMeshAndInstances = false;
    for (size_t i = 0; i < root.instances.size(); i++)
    {
        const uint32_t childSceneIndex = root.instances[i].childSceneIndex;
        if (childSceneIndex >= dag.scenes.size())
        {
            std::fprintf(stderr, "NestedInstancing: child scene index out of range in root\n");
            return false;
        }

        const ybi::USDBuildScene &child = dag.scenes[childSceneIndex];
        if (!child.meshes.empty() && !child.instances.empty())
        {
            foundPrototypeWithMeshAndInstances = true;
            break;
        }
    }
    if (!foundPrototypeWithMeshAndInstances)
    {
        std::fprintf(stderr,
                     "NestedInstancing: expected at least one referenced prototype scene with "
                     "both meshes and instances\n");
        return false;
    }

    std::unordered_map<std::string, size_t> sceneIndexByPath;
    sceneIndexByPath.reserve(dag.scenes.size());
    for (size_t i = 0; i < dag.scenes.size(); i++)
    {
        sceneIndexByPath.emplace(dag.scenes[i].path, i);
    }

    auto withInstanceIt = sceneIndexByPath.find("/PIProtoWithInstance");
    if (withInstanceIt == sceneIndexByPath.end())
    {
        std::fprintf(stderr, "NestedInstancing: /PIProtoWithInstance scene missing in DAG\n");
        return false;
    }
    const ybi::USDBuildScene &withInstanceScene = dag.scenes[withInstanceIt->second];
    if (withInstanceScene.instances.empty())
    {
        std::fprintf(stderr, "NestedInstancing: expected /PIProtoWithInstance to contain instances\n");
        return false;
    }

    auto plainIt = sceneIndexByPath.find("/PIProtoPlain");
    if (plainIt == sceneIndexByPath.end())
    {
        std::fprintf(stderr, "NestedInstancing: /PIProtoPlain scene missing in DAG\n");
        return false;
    }

    return true;
}

} // namespace

int main()
{
    int failed = 0;

    if (!RunNoInstancesTest())
    {
        failed++;
    }
    if (!RunVisibilityCullTest())
    {
        failed++;
    }
    if (!RunNormalInstancingTest())
    {
        failed++;
    }
    if (!RunNestedInstancingTest())
    {
        failed++;
    }
    if (!RunPointInstancerTransformLoadTest())
    {
        failed++;
    }

    if (failed != 0)
    {
        std::fprintf(stderr, "usd_instance_dag_test: %d test(s) failed\n", failed);
        return 1;
    }

    std::printf("usd_instance_dag_test: all tests passed\n");
    return 0;
}
