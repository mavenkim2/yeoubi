#include "scene/device_scene.h"
#include <memory>

namespace ybi
{

class DeviceArena;

struct DataContainer;
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

struct RenderSystem
{
    Scene *rootScene;
    DataContainer *dataContainer;

    DeviceScene deviceScene;

    Device *device;
    std::unique_ptr<DeviceArena> deviceArena;

    // Data managers
    GeometryManager geometryManager;

    RenderSystem(Device *device, DataContainer *container, Scene *rootScene);
    void Update();
};

} // namespace ybi
