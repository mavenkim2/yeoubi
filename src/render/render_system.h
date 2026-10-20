#include "core/arena.h"
#include "device/device_memory.h"
#include "device/device_texture.h"
#include "render/virtual_texture/manager.h"
#include "kernel/scene.h"
#include "material/material.h"
#include "material/material_compilation.h"
#include <memory>
#include <string>
#include <unordered_map>
#include <vector>

namespace ybi
{

class DeviceArena;
class VirtualTextureManager;

struct BufferHandle;
struct DataContainer;
struct OpenPbrMaterial;
struct SceneImportContext;
struct DisplayDevice;
struct Device;
struct Scene;

struct GeometryManager
{
    void BuildClusterMeshBVHsAndUpload(Device *device,
                                       DeviceArena *deviceArena,
                                       DataContainer *container,
                                       DeviceScene *deviceScene,
                                       const Array<BufferHandle> &meshHandles,
                                       Array<Object> &objects,
                                       Array<int> &objectBVHIndices);
    void BuildTriangleMeshBVHsAndUpload(Device *device,
                                        DeviceArena *deviceArena,
                                        DataContainer *container,
                                        DeviceScene *deviceScene,
                                        const Array<BufferHandle> &meshHandles,
                                        Array<Object> &objects,
                                        Array<int> &objectBVHIndices);
    void Update(Device *device,
                DeviceArena *deviceArena,
                DataContainer *container,
                Scene *rootScene,
                DeviceScene *deviceScene);
};

struct AttributeManager
{
    AttributeId counter;
    std::unordered_map<std::string, AttributeId> attributeIdMap;

    AttributeId GetOrAdd(const std::string &name);
    void FindPrimvars(DataContainer *container, SceneImportContext *context);
};

struct MaterialManager
{
    Arena arena;

    // Maps material node to MaterialTextureHandle
    std::unordered_map<const MaterialNode *, MaterialTextureHandle> textureNodeToHandle;
    VirtualTextureManager *virtualTextureManager = nullptr;
    AttributeManager attributeManager;
    MaterialCompiler materialCompiler;

    ~MaterialManager();
    void AddTextures(MaterialGraph &graph, MaterialNodeId nodeId);
    DeviceMaterialExpression ProcessMaterialGraph(OpenPbrMaterial &material,
                                                  MaterialNodeId member,
                                                  std::vector<uint32_t> &nodes);
    std::string ProcessMaterialXMaterial(MaterialGraph &material);
    void Update(Device *device,
                DisplayDevice *displayDevice,
                DataContainer *container,
                DeviceScene *deviceScene,
                SceneImportContext *importContext);
};

struct RenderSystem
{
    struct EnvmapLightResources
    {
        DeviceTexture texture;
        DeviceBuffer<float> conditionalCDF;
        DeviceBuffer<float> marginalCDF;
    };

    Scene *rootScene;
    DataContainer *dataContainer;

    DeviceScene deviceScene{};

    Device *device;
    DisplayDevice *displayDevice;
    std::unique_ptr<DeviceArena> deviceArena;

    // Data managers
    GeometryManager geometryManager;
    MaterialManager materialManager;
    std::vector<EnvmapLightResources> envmapLightResources;

    RenderSystem(Device *device,
                 DisplayDevice *displayDevice,
                 DataContainer *container,
                 Scene *rootScene,
                 VirtualTextureManager *virtualTextureManager);
    ~RenderSystem();
    void Update(SceneImportContext *context);
    void UpdateEnvmapLights();
    void ReleaseEnvmapLights();
};

} // namespace ybi
