#include "render/scene_memory.h"

#include <cstdint>

namespace ybi
{
namespace
{

template <typename T>
uint64_t ArrayBytes(const Array<T> &array)
{
    return static_cast<uint64_t>(array.size()) * sizeof(T);
}

uint64_t StringBytes(const std::string &value)
{
    return static_cast<uint64_t>(value.size());
}

uint64_t AttributeBytes(const Attribute &attribute)
{
    return StringBytes(attribute.name) + ArrayBytes(attribute.data) + ArrayBytes(attribute.indices);
}

uint64_t MaterialTextureInputBytes(const MaterialTextureInput &input)
{
    return StringBytes(input.inputName) + StringBytes(input.texturePath) + StringBytes(input.swizzle);
}

uint64_t MaterialInfoBytes(const MaterialInfo &material)
{
    uint64_t bytes = StringBytes(material.materialPath) + StringBytes(material.ntcDiffuseFile) +
                     StringBytes(material.ntcDiffuseTextureName) + sizeof(material.packed);
    for (const MaterialTextureInput &input : material.textureInputs)
    {
        bytes += MaterialTextureInputBytes(input);
    }
    return bytes;
}

uint64_t LightInfoBytes(const LightInfo &light)
{
    return StringBytes(light.lightPath) + StringBytes(light.texturePath) + sizeof(light.packed);
}

uint64_t CurvesHostBytes(const Curves &curves)
{
    return ArrayBytes(curves.GetVertices()) + ArrayBytes(curves.GetWidths()) +
           static_cast<uint64_t>(curves.GetNumCurves() + 1u) * sizeof(int);
}

uint64_t MicropolygonMeshBytes(const MicropolygonMesh &mesh)
{
    return ArrayBytes(mesh.positions) + ArrayBytes(mesh.grids);
}

uint64_t SubdivisionMeshBytes(const SubdivisionMesh &mesh)
{
    uint64_t bytes = StringBytes(mesh.primPath) + ArrayBytes(mesh.vertices) + ArrayBytes(mesh.indices) +
                     ArrayBytes(mesh.vertsPerFace) + ArrayBytes(mesh.cornerIndices) +
                     ArrayBytes(mesh.cornerSharpnesses) + ArrayBytes(mesh.creaseIndices) +
                     ArrayBytes(mesh.creaseLengths) + ArrayBytes(mesh.creaseSharpnesses) +
                     ArrayBytes(mesh.holeIndices);
    for (const Attribute &attribute : mesh.attributes)
    {
        bytes += AttributeBytes(attribute);
    }
    return bytes;
}

} // namespace

uint64_t ComputeMeshHostBytes(const Mesh &mesh)
{
    uint64_t bytes = ArrayBytes(mesh.positions) + ArrayBytes(mesh.indices);
    for (const Attribute &attribute : mesh.attributes)
    {
        bytes += AttributeBytes(attribute);
    }
    return bytes;
}

uint64_t ComputeMeshVectorHostBytes(const std::vector<Mesh> &meshes)
{
    uint64_t bytes = 0u;
    for (const Mesh &mesh : meshes)
    {
        bytes += ComputeMeshHostBytes(mesh);
    }
    return bytes;
}

uint64_t ComputeSceneHostBytes(const Scene &scene)
{
    uint64_t bytes = ArrayBytes(scene.primitives) + ArrayBytes(scene.primitiveCollections) +
                     ArrayBytes(scene.attributes) +
                     static_cast<uint64_t>(scene.instances.size()) * sizeof(Instance) +
                     static_cast<uint64_t>(scene.childScenes.size()) * sizeof(Scene *);
    for (const Mesh &mesh : scene.meshes)
    {
        bytes += ComputeMeshHostBytes(mesh);
    }
    for (const Curves &curves : scene.curves)
    {
        bytes += CurvesHostBytes(curves);
    }
    for (const MicropolygonMesh &mesh : scene.micropolygonMeshes)
    {
        bytes += MicropolygonMeshBytes(mesh);
    }
    for (const SubdivisionMesh &mesh : scene.subdivisionMeshes)
    {
        bytes += SubdivisionMeshBytes(mesh);
    }
    for (const Attribute &attribute : scene.attributes)
    {
        bytes += AttributeBytes(attribute);
    }
    return bytes;
}

uint64_t ComputeScenePoolHostBytes(const ScenePool &scenePool)
{
    uint64_t bytes = StringBytes(scenePool.camera.path);
    for (const std::unique_ptr<Scene> &scene : scenePool.scenes)
    {
        if (scene)
        {
            bytes += ComputeSceneHostBytes(*scene);
        }
    }
    for (const MaterialInfo &material : scenePool.materials)
    {
        bytes += MaterialInfoBytes(material);
    }
    for (const LightInfo &light : scenePool.lights)
    {
        bytes += LightInfoBytes(light);
    }
    return bytes;
}

} // namespace ybi
