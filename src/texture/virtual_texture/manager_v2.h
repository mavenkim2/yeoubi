#include "device/device.h"
#include "texture/virtual_texture/manager.h"
#include "util/array.h"

namespace ybi
{

struct FeedbackV2
{
    uint32_t udimID;
};

struct VirtualTextureManagerV2
{
    struct TextureState
    {
        std::string tileFilePath;
    };

    struct UDIMData
    {
        uint32_t udim;
        uint32_t textureStateID;
    };

    // Hash map again
    // this one will be different
    // this will be interactive vis i guess

    Device *device;
    // Device textures
    DeviceTexture texturePool[2];

    // Maps from one per texture ID to one per UDIM ID
    Array<uint32_t> textureUdimIDOffsets;
    Array<UDIMData> textureUdimIDData;

    Array<TextureState> textureStates;

    // Feedback
    DeviceMemoryView<uint32_t> udimBitArray;
    uint32_t *hostUdimBitArray;

    bool RegisterTexture(const texture::VirtualTextureRegisterInput &input);
    void ProcessFeedback();
};

bool VirtualTextureManagerV2::RegisterTexture(const texture::VirtualTextureRegisterInput &input)
{
    const uint32_t textureID = input.textureId;

    TextureState state;
    state.tileFilePath = input.tileFilePath;

    if (textureID >= textureStates.size())
    {
        textureStates.Resize(textureID + 1);
    }
    textureStates[textureID] = state;

    uint32_t currentOffset = uint32_t(textureUdimIDData.size());
    textureUdimIDOffsets.EmplaceBack(currentOffset);

    textureUdimIDData.Resize(textureUdimIDData.size() + input.activeUdims.size());

    for (uint32_t index = currentOffset; index < uint32_t(textureUdimIDData.size()); index++)
    {
        UDIMData data;
        data.udim = input.activeUdims[index - currentOffset];
        data.textureStateID = textureID;

        textureUdimIDData[index] = data;
    }
}

inline uint32_t FirstBitLow(uint32_t x)
{
#error something
    return 0;
}

void VirtualTextureManagerV2::ProcessFeedback()
{
    const uint32_t numEntries = uint32_t(udimBitArray.size());

    if (device->GetKind() != DeviceKind::CPU)
    {
        device->CopyBytesToHost(hostUdimBitArray, udimBitArray, udimBitArray.numBytes());
    }

    uint32_t *bitArray =
        device->GetKind() == DeviceKind::CPU ? udimBitArray.data() : hostUdimBitArray;

    for (uint32_t entryIndex = 0; entryIndex < numEntries; entryIndex++)
    {
        uint32_t entry = bitArray[entryIndex];
        while (entry)
        {
            uint32_t bitIndex = FirstBitLow(entry);
            uint32_t descriptorIndex = 32 * entryIndex + bitIndex;
            entry &= entry - 1;
        }
    }
}

} // namespace ybi
