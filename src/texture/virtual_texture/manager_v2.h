#include "device/device.h"
#include "texture/texture_format.h"
#include "texture/virtual_texture/manager.h"
#include "texture/virtual_texture/tile_file.h"
#include "util/array.h"
#include "util/assert.h"
#include "util/platform.h"
#include <mutex>

namespace ybi
{

struct FeedbackV2
{
    uint32_t udimID;
    uint32_t mip;
    uint32_t tileX;
    uint32_t tileY;
};

template <typename T>
struct RingBuffer
{
    Array<T> entries;
    uint64_t readOffset;
    uint64_t writeOffset;

    uint32_t max;

    RingBuffer(Arena *arena, uint32_t max);
    bool Write(T *vals, uint32_t num);
    void WriteWithOverwrite(T *vals, uint32_t num);
    void SynchronizedWrite(std::mutex *mutex, T *vals, uint32_t num);
    T *Read(Arena *arena, uint32_t &num);
    T *SynchronizedRead(std::mutex *mutex, Arena *arena, uint32_t &num);
};

struct VirtualTextureManagerV2
{
    struct TextureState
    {
        texture::VirtualTextureTileFile tileFile;
        std::string tileFilePath;
        TextureFormat format;
    };

    struct UDIMData
    {
        uint32_t udim;
        uint32_t textureStateID;
    };

    struct TexturePool
    {
        Array<DeviceTextureV2> textures;
        Array<std::pair<uint32_t, uint32_t>> freeList;
    };

    // Hash map again
    // this one will be different
    // this will be interactive vis i guess

    Device *device;
    // Device textures
    std::array<TexturePool, uint32_t(TextureFormat::MAX)> texturePools;

    RingBuffer ringBuffer;

    // Maps from one per texture ID to one per UDIM ID
    Array<uint32_t> textureUdimIDOffsets;
    Array<UDIMData> textureUdimIDData;

    Array<TextureState> textureStates;

    // Feedback
    DeviceMemoryView<uint32_t> udimBitArray;
    uint32_t *hostUdimBitArray;

    DeviceMemoryView<FeedbackV2> feedbackArray;
    FeedbackV2 *hostFeedbackArray;

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

void VirtualTextureManagerV2::ProcessFeedback()
{
    const uint32_t numEntries = uint32_t(udimBitArray.size());

    device->CopyBytesToHost(hostUdimBitArray, udimBitArray, udimBitArray.numBytes());
    device->CopyBytesToHost(hostFeedbackArray, feedbackArray, feedbackArray.numBytes());

    uint32_t *bitArray = hostUdimBitArray;
    FeedbackV2 *feedback = hostFeedbackArray;

    // TODO: compact feedback?
    uint32_t count = feedback[0].udimID;
    for (uint32_t feedbackIndex = 1; feedbackIndex < count + 1; feedbackIndex++)
    {
        // TODO: something
        const FeedbackV2 &f = feedback[feedbackIndex];
        const UDIMData &data = textureUdimIDData[f.udimID];

        const TextureState &textureState = textureStates[data.textureStateID];
        const std::string &tileFilePath = textureState.tileFilePath;

        TextureFormat format;
        uint32_t width, height;
        uint64_t outSourceBytes;
        std::string error;
        std::vector<unsigned char> outPixels;
        texture::ReadVirtualTextureTile(&textureState.tileFilePath,
                                        data.udim,
                                        f.mip,
                                        f.tileX,
                                        f.tileY,
                                        &outPixels,
                                        &format,
                                        &width,
                                        &height,
                                        &outSourceBytes,
                                        &error);
        YBI_ASSERT(format == textureState.format);

        TexturePool &pool = texturePools[uint32_t(format)];
        if (pool.freeList.size())
        {
            std::pair<uint32_t, uint32_t> poolKey = pool.freeList.Pop();
            DeviceTextureV2 &texture = pool.textures[poolKey.first];
        }
    }

#if 0
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
#endif
}

} // namespace ybi
