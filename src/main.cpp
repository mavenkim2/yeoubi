#include "core/arena.h"
#include "core/assert.h"
#include "device/device.h"
#include "device/kernel_params.h"
#include "io/usd.h"
#include "scene/scene.h"
#include <memory>
#include <string>

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include <stb_image_write.h>

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

    YBI_ASSERT(scene->container->GetBVHs(GeometryType::Scene).size() == 1);

    AffineTransform renderFromCamera;
    Matrix4x4 cameraFromRaster;

    const Camera &camera = scene->container->camera;
    bool result = camera.cameraFromWorld.TryInverse(&renderFromCamera);
    YBI_ASSERT(result);
    result = camera.rasterFromCamera.TryInverse(&cameraFromRaster);
    YBI_ASSERT(result);

    // Image
    DeviceBuffer imageHandle = device->Alloc(camera.width * camera.height * 4);

    KernelParams params;
    params.bvhHandle = scene->container->GetBVHs(GeometryType::Scene)[0].handle.handle;
    params.imageHandle = reinterpret_cast<uint64_t>(imageHandle.Data());
    params.cameraFromRaster = cameraFromRaster;
    params.renderFromCamera = renderFromCamera;
    params.maxDepth = 3;
    params.lensRadius = 0.f; // TODO

    DeviceBuffer cudaParams = device->Alloc(sizeof(KernelParams));
    device->CopyToDevice(cudaParams.Data(), &params, sizeof(KernelParams));
    device->DispatchKernel(1920, 1080, cudaParams);

    void *outBytes = malloc(1920 * 1080 * 4);
    device->CopyToHost(outBytes, imageHandle.Data(), 1920 * 1080 * 4);

    stbi_write_png("test.png", 1920, 1080, 4, outBytes, 1920 * 4);

    return 0;
}
