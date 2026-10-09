#include "render/path_guiding/lut.h"

#include "core/executable_path.h"
#include "device/device.h"
#include <array>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <memory>

namespace ybi
{

VMFSplitTable LoadVMFSplitTable(Device &device)
{
    constexpr size_t byteCount = (VMFSplitTableSize + 1) * sizeof(float);
    const std::filesystem::path executableDirectory = ExecutableDirectory();
    if (executableDirectory.empty())
    {
        std::fprintf(stderr, "Unable to locate the executable directory.\n");
        return {};
    }
    const std::filesystem::path tablePath = executableDirectory / VMFSplitTableFilename;
    std::ifstream input(tablePath, std::ios::binary | std::ios::ate);
    if (!input || input.tellg() != static_cast<std::streamoff>(byteCount))
    {
        std::fprintf(
            stderr, "Failed to load the vMF split table: %s\n", tablePath.string().c_str());
        return {};
    }

    std::array<float, VMFSplitTableSize + 1> values;
    input.seekg(0, std::ios::beg);
    input.read(reinterpret_cast<char *>(values.data()), byteCount);
    if (!input)
    {
        std::fprintf(
            stderr, "Failed to read the vMF split table: %s\n", tablePath.string().c_str());
        return {};
    }

    DeviceTextureInfo textureInfo{};
    textureInfo.format = TextureFormat::R32_FLOAT;
    textureInfo.filterType = FilterType::Linear;
    textureInfo.width = VMFSplitTableSize;
    VMFSplitTable table;
    table.maxInvSqrtK = values[0];
    table.texture = device.CreateTexture(textureInfo);

    std::unique_ptr<DeviceStream> uploadStream = device.CreateDeviceStream();
    DeviceTextureCopyDesc copy{};
    copy.src = values.data() + 1;
    copy.srcRowPitch = VMFSplitTableSize * sizeof(float);
    copy.srcSlicePitch = copy.srcRowPitch;
    copy.widthInBytes = copy.srcRowPitch;
    copy.height = 1;
    copy.depth = 1;
    device.CopyToDevice3DAsync(*uploadStream, table.texture, copy);
    device.SyncDeviceStream(*uploadStream);
    return table;
}

} // namespace ybi
