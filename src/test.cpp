#include "device/cuda_device.h"
#include "device/device.h"
#include "io/usd/load.h"
#include "scene/scene.h"
#include "util/float3.h"
#include <cassert>
#include <chrono>
#include <cstdlib>
#include <cstdio>
#include <memory>
#include <string>
#include <vector>

using namespace ybi;

namespace
{
struct CliOptions
{
    DeviceKind deviceKind = DeviceKind::GPU;
    std::string usdFilePath;
};

static void PrintUsage(const char *exe)
{
    printf("usage: %s <-gpu|-cpu> <usd_file>\n", exe);
}

static CliOptions ParseCli(int argc, char **argv)
{
    CliOptions options = {};
    for (int i = 1; i < argc; i++)
    {
        const std::string arg = argv[i];
        if (arg == "-gpu")
        {
            options.deviceKind = DeviceKind::GPU;
            continue;
        }
        if (arg == "-cpu")
        {
            options.deviceKind = DeviceKind::CPU;
            continue;
        }
        if (!arg.empty() && arg[0] == '-')
        {
            PrintUsage(argv[0]);
            std::exit(1);
        }
        if (!options.usdFilePath.empty())
        {
            PrintUsage(argv[0]);
            std::exit(1);
        }
        options.usdFilePath = arg;
    }

    if (options.usdFilePath.empty())
    {
        PrintUsage(argv[0]);
        std::exit(1);
    }
    return options;
}

static Device *CreateDevice(DeviceKind kind, std::unique_ptr<Device> &storage)
{
    if (kind == DeviceKind::CPU)
    {
        fprintf(stderr, "CPU device backend not implemented yet.\n");
        return nullptr;
    }

#if defined(WITH_CUDA) && defined(WITH_OPTIX)
    storage = std::make_unique<CUDADevice>();
    return storage.get();
#else
    (void)storage;
    fprintf(stderr, "GPU backend requires WITH_CUDA and WITH_OPTIX.\n");
    return nullptr;
#endif
}
} // namespace

int main(int argc, char **argv)
{
    const CliOptions options = ParseCli(argc, argv);

    ScenePool scenePool;
    LoadUSDScene(&scenePool, options.usdFilePath);
    if (scenePool.scenes.empty() || scenePool.rootSceneIndex >= scenePool.scenes.size())
    {
        fprintf(stderr, "Failed to load USD scene or invalid root: %s\n", options.usdFilePath.c_str());
        return 1;
    }

    std::unique_ptr<Device> deviceStorage;
    Device *device = CreateDevice(options.deviceKind, deviceStorage);
    if (!device)
    {
        return 1;
    }

    ScenePool flattenedScenePool;
    std::string flattenError;
    if (!FlattenScenePoolToRootChildren(&scenePool, &flattenedScenePool, &flattenError))
    {
        fprintf(stderr, "Failed to flatten USD ScenePool: %s\n", flattenError.c_str());
        return 1;
    }
    if (flattenedScenePool.scenes.empty() ||
        flattenedScenePool.rootSceneIndex >= flattenedScenePool.scenes.size())
    {
        fprintf(stderr, "Flattened ScenePool invalid: %s\n", options.usdFilePath.c_str());
        return 1;
    }

    std::vector<SceneMeshUploadRef> meshUploadRefs;
    CollectScenePoolMeshUploadRefs(&flattenedScenePool, &meshUploadRefs);

    printf("bvh start\n");
    auto start = std::chrono::high_resolution_clock::now();
    for (const std::unique_ptr<Scene> &scenePtr : flattenedScenePool.scenes)
    {
        Scene *scene = scenePtr.get();
        assert(scene);
        device->BuildBVH(scene);
    }
    auto end = std::chrono::high_resolution_clock::now();
    printf("bvh end\n");

    std::vector<DeviceMemoryView<uint8_t>> ownedMeshBuffers;
    ownedMeshBuffers.reserve(meshUploadRefs.size() * 2);
    for (const SceneMeshUploadRef &uploadRef : meshUploadRefs)
    {
        assert(uploadRef.mesh);
        const Mesh &mesh = *uploadRef.mesh;
        if (mesh.positions.size() == 0 || mesh.indices.size() == 0)
        {
            continue;
        }

        const size_t positionsBytes = sizeof(ybi::float3) * mesh.positions.size();
        const size_t indicesBytes = sizeof(int) * mesh.indices.size();

        DeviceMemoryView<ybi::float3> positionsBuffer =
            device->Alloc<ybi::float3>(mesh.positions.size());
        DeviceMemoryView<int> indicesBuffer = device->Alloc<int>(mesh.indices.size());
        device->CopyBytesToDevice(Device::ByteView(positionsBuffer),
                                  mesh.positions.data(),
                                  positionsBytes);
        device->CopyBytesToDevice(
            Device::ByteView(indicesBuffer), mesh.indices.data(), indicesBytes);

        ownedMeshBuffers.push_back(Device::ByteView(positionsBuffer));
        ownedMeshBuffers.push_back(Device::ByteView(indicesBuffer));
    }

    for (DeviceMemoryView<uint8_t> &buffer : ownedMeshBuffers)
    {
        device->FreeBytes(buffer);
    }

    const std::chrono::duration<double, std::milli> elapsed = end - start;
    printf("time elapsed: %f ms\n", elapsed.count());
    printf("bvh size: %zi\n", device->GetBVHAllocatedBytes());

    return 0;
}
