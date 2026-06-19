#include "core/arena.h"
#include "core/assert.h"
#include "device/device.h"
#include "io/usd.h"
#include "scene/scene.h"
#include <memory>
#include <string>

using namespace ybi;

int main(int argc, char **argv)
{
    YBI_ASSERT(argc == 2);

    // char *filename = argv[0];

    std::string filename(argv[1]);
    // "D:\\Distributable_2023_Davinci_v2\\Distributable_2023_Davinci\\shot\\rt\\rt_base.usda";
    // "D:\\ALab-2.2.0\\ALab\\entry.usda";
    Arena arena;
    DataContainer container;

    std::unique_ptr<SceneLoader> sceneLoader = SceneLoader::CreateSceneLoader(filename);
    YBI_ASSERT(sceneLoader);
    Scene *scene = sceneLoader->LoadScene(&arena, &container);

    std::unique_ptr<Device> device = Device::CreateDevice(DeviceKind::GPU);

    RenderSystem system(device.get(), &container, scene);
    system.Update();

    return 0;
}
