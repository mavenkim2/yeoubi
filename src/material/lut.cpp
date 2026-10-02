#include "material/lut.h"

#include "core/executable_path.h"
#include "device/device.h"
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <memory>
#include <vector>

namespace ybi
{

constexpr uint32_t OpenPbrEnergyTableSize = 32;
constexpr const char *OpenPbrEnergyTableFilename =
    "openpbr_opaque_dielectric_energy_complement.bin";
constexpr const char *OpenPbrAverageEnergyTableFilename =
    "openpbr_opaque_dielectric_average_energy_complement.bin";

DeviceTexture LoadOpenPbrOpaqueDielectricEnergyComplementTexture(Device &device)
{
    constexpr size_t valueCount =
        OpenPbrEnergyTableSize * OpenPbrEnergyTableSize * OpenPbrEnergyTableSize;
    constexpr size_t byteCount = valueCount * sizeof(uint16_t);

    const std::filesystem::path executableDirectory = ExecutableDirectory();
    if (executableDirectory.empty())
    {
        std::fprintf(stderr, "Unable to locate the executable directory.\n");
        return {};
    }
    const std::filesystem::path tablePath = executableDirectory / OpenPbrEnergyTableFilename;
    std::ifstream input(tablePath, std::ios::binary | std::ios::ate);
    if (!input || input.tellg() != static_cast<std::streamoff>(byteCount))
    {
        std::fprintf(stderr,
                     "Failed to load the 32x32x32 OpenPBR energy table: %s\n",
                     tablePath.string().c_str());
        return {};
    }

    std::vector<uint16_t> values(valueCount);
    input.seekg(0, std::ios::beg);
    input.read(reinterpret_cast<char *>(values.data()), byteCount);
    if (!input)
    {
        std::fprintf(
            stderr, "Failed to read the OpenPBR energy table: %s\n", tablePath.string().c_str());
        return {};
    }

    DeviceTextureInfo textureInfo{};
    textureInfo.format = TextureFormat::R16_UNORM;
    textureInfo.filterType = FilterType::Linear;
    textureInfo.width = OpenPbrEnergyTableSize;
    textureInfo.height = OpenPbrEnergyTableSize;
    textureInfo.depth = OpenPbrEnergyTableSize;
    DeviceTexture texture = device.CreateTexture(textureInfo);

    std::unique_ptr<DeviceStream> uploadStream = device.CreateDeviceStream();
    DeviceTextureCopyDesc copy{};
    copy.src = values.data();
    copy.srcRowPitch = OpenPbrEnergyTableSize * sizeof(uint16_t);
    copy.srcSlicePitch = OpenPbrEnergyTableSize * copy.srcRowPitch;
    copy.widthInBytes = copy.srcRowPitch;
    copy.height = OpenPbrEnergyTableSize;
    copy.depth = OpenPbrEnergyTableSize;
    device.CopyToDevice3DAsync(*uploadStream, texture, copy);
    device.SyncDeviceStream(*uploadStream);
    return texture;
}

DeviceTexture LoadOpenPbrOpaqueDielectricAverageEnergyComplementTexture(Device &device)
{
    constexpr size_t valueCount = OpenPbrEnergyTableSize * OpenPbrEnergyTableSize;
    constexpr size_t byteCount = valueCount * sizeof(uint16_t);

    const std::filesystem::path executableDirectory = ExecutableDirectory();
    if (executableDirectory.empty())
    {
        std::fprintf(stderr, "Unable to locate the executable directory.\n");
        return {};
    }
    const std::filesystem::path tablePath =
        executableDirectory / OpenPbrAverageEnergyTableFilename;
    std::ifstream input(tablePath, std::ios::binary | std::ios::ate);
    if (!input || input.tellg() != static_cast<std::streamoff>(byteCount))
    {
        std::fprintf(stderr,
                     "Failed to load the 32x32 OpenPBR average energy table: %s\n",
                     tablePath.string().c_str());
        return {};
    }

    std::vector<uint16_t> values(valueCount);
    input.seekg(0, std::ios::beg);
    input.read(reinterpret_cast<char *>(values.data()), byteCount);
    if (!input)
    {
        std::fprintf(stderr,
                     "Failed to read the OpenPBR average energy table: %s\n",
                     tablePath.string().c_str());
        return {};
    }

    DeviceTextureInfo textureInfo{};
    textureInfo.format = TextureFormat::R16_UNORM;
    textureInfo.filterType = FilterType::Linear;
    textureInfo.width = OpenPbrEnergyTableSize;
    textureInfo.height = OpenPbrEnergyTableSize;
    textureInfo.layers = 1;
    DeviceTexture texture = device.CreateTexture(textureInfo);

    std::unique_ptr<DeviceStream> uploadStream = device.CreateDeviceStream();
    DeviceTextureCopyDesc copy{};
    copy.src = values.data();
    copy.srcRowPitch = OpenPbrEnergyTableSize * sizeof(uint16_t);
    copy.srcSlicePitch = OpenPbrEnergyTableSize * copy.srcRowPitch;
    copy.widthInBytes = copy.srcRowPitch;
    copy.height = OpenPbrEnergyTableSize;
    copy.depth = 1;
    device.CopyToDevice3DAsync(*uploadStream, texture, copy);
    device.SyncDeviceStream(*uploadStream);
    return texture;
}

} // namespace ybi
