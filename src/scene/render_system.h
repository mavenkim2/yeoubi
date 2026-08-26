#include "core/arena.h"
#include "scene/device_scene.h"
#include "scene/material.h"
#include <memory>

namespace ybi
{

class DeviceArena;
class VirtualTextureManager;

struct DataContainer;
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

struct MaterialManager
{
    Arena arena;
    std::unique_ptr<VirtualTextureManager> virtualTextureManager;

    ~MaterialManager();
    void AddTextures(OpenPbrMaterial &material, MaterialNodeId member);
    void ProcessMaterialGraph(OpenPbrMaterial &material,
                              const char *memberName,
                              MaterialNodeId member,
                              std::vector<uint32_t> &nodes);
    void Update(Device *device,
                DisplayDevice *displayDevice,
                DataContainer *container,
                DeviceScene *deviceScene);
};

struct RenderSystem
{
    Scene *rootScene;
    DataContainer *dataContainer;

    DeviceScene deviceScene;

    Device *device;
    std::unique_ptr<DeviceArena> deviceArena;

    // Data managers
    GeometryManager geometryManager;
    MaterialManager materialManager;

    RenderSystem(Device *device, DataContainer *container, Scene *rootScene);
    ~RenderSystem();
    void Update();
};

} // namespace ybi
