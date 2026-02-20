#include "io/usd/instance_dag_build.h"

#include <pxr/usd/sdf/layer.h>
#include <pxr/usd/usd/stage.h>

#include <cstdio>
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
    if (!RunNormalInstancingTest())
    {
        failed++;
    }
    if (!RunNestedInstancingTest())
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
