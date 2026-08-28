#include "core/arena.h"
#include "render/virtual_texture_manager.h"
#include "scene/device_scene.h"
#include "scene/material.h"
#include <memory>
#include <unordered_map>

namespace ybi
{

class DeviceArena;
class VirtualTextureManager;

struct DataContainer;
struct SceneImportContext;
struct DisplayDevice;
struct Device;
struct Scene;

struct GeometryManager
{
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
    std::unique_ptr<VirtualTextureManager> virtualTextureManager;
    AttributeManager attributeManager;

    ~MaterialManager();
    void AddTextures(OpenPbrMaterial &material, MaterialNodeId member);
    void ProcessMaterialGraph(OpenPbrMaterial &material,
                              const char *memberName,
                              MaterialNodeId member,
                              std::vector<uint32_t> &nodes);
    void Update(Device *device,
                DisplayDevice *displayDevice,
                DataContainer *container,
                DeviceScene *deviceScene,
                SceneImportContext *importContext);
};

struct RenderSystem
{
    Scene *rootScene;
    DataContainer *dataContainer;

    DeviceScene deviceScene;

    Device *device;
    DisplayDevice *displayDevice;
    std::unique_ptr<DeviceArena> deviceArena;

    // Data managers
    GeometryManager geometryManager;
    MaterialManager materialManager;

    RenderSystem(Device *device,
                 DisplayDevice *displayDevice,
                 DataContainer *container,
                 Scene *rootScene);
    ~RenderSystem();
    void Update(SceneImportContext *context);
};

} // namespace ybi
