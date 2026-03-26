#include "texture/virtual_texture/manager.h"

#include <algorithm>
#include <cmath>

namespace ybi
{
namespace texture
{

namespace
{
static constexpr uint32_t kStreamPagesPerTexture = 256u;

uint64_t ComputeTileFileHostBytes(const VirtualTextureTileFile &file)
{
    uint64_t bytes = static_cast<uint64_t>(file.path.size());
    for (const auto &udimEntry : file.udims)
    {
        bytes += sizeof(uint32_t) + sizeof(VirtualTextureUdimTable);
        for (const VirtualTextureMipTable &mip : udimEntry.second.mips)
        {
            bytes += sizeof(VirtualTextureMipTable);
            bytes += static_cast<uint64_t>(mip.records.size()) * sizeof(VirtualTextureTileRecord);
        }
    }
    return bytes;
}

} // namespace

bool VirtualTextureManager::AllocateDeviceState(std::string *outError)
{
    YBI_ASSERT(device_);
    if (!pageTableEntriesHost_.empty())
    {
        const size_t bytes = pageTableEntriesHost_.size() * sizeof(uint32_t);
        pageTableEntriesDevice_ = device_->AllocBytes(bytes);
        device_->CopyBytesToDevice(pageTableEntriesDevice_, pageTableEntriesHost_.data(), bytes);
    }
    if (!pageTableMipOffsets_.empty())
    {
        const size_t bytes = pageTableMipOffsets_.size() * sizeof(uint32_t);
        pageTableMipOffsetsDevice_ = device_->AllocBytes(bytes);
        device_->CopyBytesToDevice(pageTableMipOffsetsDevice_, pageTableMipOffsets_.data(), bytes);
    }
    if (!pageTableMipWidths_.empty())
    {
        const size_t bytes = pageTableMipWidths_.size() * sizeof(uint32_t);
        pageTableMipWidthsDevice_ = device_->AllocBytes(bytes);
        device_->CopyBytesToDevice(pageTableMipWidthsDevice_, pageTableMipWidths_.data(), bytes);
    }
    if (!pageTableMipHeights_.empty())
    {
        const size_t bytes = pageTableMipHeights_.size() * sizeof(uint32_t);
        pageTableMipHeightsDevice_ = device_->AllocBytes(bytes);
        device_->CopyBytesToDevice(pageTableMipHeightsDevice_, pageTableMipHeights_.data(), bytes);
    }

    if (!mipInfosHost_.empty())
    {
        const size_t bytes = mipInfosHost_.size() * sizeof(LaunchParams::VirtualTextureMipInfo);
        mipInfosDevice_ = device_->AllocBytes(bytes);
        device_->CopyBytesToDevice(mipInfosDevice_, mipInfosHost_.data(), bytes);
    }
    if (!udimInfosHost_.empty())
    {
        const size_t bytes = udimInfosHost_.size() * sizeof(LaunchParams::VirtualTextureUdimInfo);
        udimInfosDevice_ = device_->AllocBytes(bytes);
        device_->CopyBytesToDevice(udimInfosDevice_, udimInfosHost_.data(), bytes);
    }

    for (TextureState &texture : textures_)
    {
        if (texture.tailPixelsHost.empty())
        {
            continue;
        }
        texture.tailPixelsDevice = device_->AllocBytes(texture.tailPixelsHost.size());
        device_->CopyBytesToDevice(
            texture.tailPixelsDevice, texture.tailPixelsHost.data(), texture.tailPixelsHost.size());
        textureMetaHost_[texture.textureId].tailPixels =
            reinterpret_cast<unsigned long long>(texture.tailPixelsDevice.data());
    }

    if (!textureMetaHost_.empty())
    {
        const size_t bytes = textureMetaHost_.size() * sizeof(LaunchParams::VirtualTextureTextureMeta);
        textureMetaDevice_ = device_->AllocBytes(bytes);
        device_->CopyBytesToDevice(textureMetaDevice_, textureMetaHost_.data(), bytes);
    }

    const size_t pageBytes =
        static_cast<size_t>(config_.pageSize) * static_cast<size_t>(config_.pageSize) * 4u;
    streamPageCountX_ = kStreamPagesPerTexture;
    streamPageCountY_ = 0u;
    streamPageCount_ = 0u;
    streamTextureBytes_ = static_cast<size_t>(streamPageCountX_) * pageBytes;
    maxStreamTextureCount_ =
        static_cast<uint32_t>(config_.cacheBytes / std::max<size_t>(1u, streamTextureBytes_));

    streamTextures_.clear();
    streamTextureHandlesHost_.assign(maxStreamTextureCount_, 0ull);
    if (!streamTextureHandlesHost_.empty())
    {
        const size_t bytes = streamTextureHandlesHost_.size() * sizeof(unsigned long long);
        streamTextureHandlesDevice_ = device_->AllocBytes(bytes);
        device_->CopyBytesToDevice(
            streamTextureHandlesDevice_, streamTextureHandlesHost_.data(), bytes);
    }
    streamSlots_.clear();
    keyToStreamSlot_.clear();
    lruSlots_.clear();
    freeSlots_.clear();
    freeSlots_.reserve(static_cast<size_t>(maxStreamTextureCount_) * static_cast<size_t>(streamPageCountX_));
    (void)outError;
    return true;
}

VirtualTextureMemoryStats VirtualTextureManager::GetMemoryStats() const
{
    VirtualTextureMemoryStats stats = {};
    stats.hostPageTableBytes += static_cast<uint64_t>(pageTableEntriesHost_.size()) * sizeof(uint32_t);
    stats.hostPageTableBytes += static_cast<uint64_t>(pageTableMipOffsets_.size()) * sizeof(uint32_t);
    stats.hostPageTableBytes += static_cast<uint64_t>(pageTableMipWidths_.size()) * sizeof(uint32_t);
    stats.hostPageTableBytes += static_cast<uint64_t>(pageTableMipHeights_.size()) * sizeof(uint32_t);

    stats.hostMetaBytes += static_cast<uint64_t>(mipInfosHost_.size()) *
                           sizeof(LaunchParams::VirtualTextureMipInfo);
    stats.hostMetaBytes += static_cast<uint64_t>(udimInfosHost_.size()) *
                           sizeof(LaunchParams::VirtualTextureUdimInfo);
    stats.hostMetaBytes += static_cast<uint64_t>(textureMetaHost_.size()) *
                           sizeof(LaunchParams::VirtualTextureTextureMeta);
    stats.devicePageTableBytes += pageTableEntriesDevice_.numBytes();
    stats.devicePageTableBytes += pageTableMipOffsetsDevice_.numBytes();
    stats.devicePageTableBytes += pageTableMipWidthsDevice_.numBytes();
    stats.devicePageTableBytes += pageTableMipHeightsDevice_.numBytes();
    stats.deviceMetaBytes += mipInfosDevice_.numBytes();
    stats.deviceMetaBytes += udimInfosDevice_.numBytes();
    stats.deviceMetaBytes += textureMetaDevice_.numBytes();
    stats.hostMetaBytes +=
        static_cast<uint64_t>(streamTextureHandlesHost_.size()) * sizeof(unsigned long long);
    stats.deviceMetaBytes += streamTextureHandlesDevice_.numBytes();
    stats.deviceStreamBytes +=
        static_cast<uint64_t>(streamTextures_.size()) * static_cast<uint64_t>(streamTextureBytes_);

    for (const TextureState &texture : textures_)
    {
        stats.hostTailBytes += static_cast<uint64_t>(texture.tailPixelsHost.size());
        if (texture.tileFileOpen)
        {
            stats.hostMetaBytes += ComputeTileFileHostBytes(texture.tileFile);
        }
        else
        {
            stats.hostMetaBytes += static_cast<uint64_t>(texture.tileFilePath.size());
        }
        stats.deviceTailBytes += texture.tailPixelsDevice.numBytes();
    }
    return stats;
}

bool VirtualTextureManager::Finalize(std::string *outError)
{
    if (!device_)
    {
        if (outError)
        {
            *outError = "VirtualTextureManager Finalize: manager not initialized";
        }
        return false;
    }
    if (finalized_)
    {
        return true;
    }
    if (!BuildMipReservations(outError))
    {
        return false;
    }
    if (!AllocateDeviceState(outError))
    {
        return false;
    }
    finalized_ = true;
    return true;
}

void VirtualTextureManager::BindLaunchParams(LaunchParams *params) const
{
    if (!params)
    {
        return;
    }
    params->virtualTextureEnabled = finalized_ ? 1 : 0;
    params->virtualTexturePageTableEntries =
        finalized_ ? reinterpret_cast<unsigned long long>(pageTableEntriesDevice_.data()) : 0ull;
    params->virtualTexturePageTableMipOffsets =
        finalized_ ? reinterpret_cast<unsigned long long>(pageTableMipOffsetsDevice_.data()) : 0ull;
    params->virtualTexturePageTableMipWidths =
        finalized_ ? reinterpret_cast<unsigned long long>(pageTableMipWidthsDevice_.data()) : 0ull;
    params->virtualTexturePageTableMipHeights =
        finalized_ ? reinterpret_cast<unsigned long long>(pageTableMipHeightsDevice_.data()) : 0ull;
    params->virtualTexturePageTableMipCount = finalized_ ? static_cast<int>(pageTableMipOffsets_.size()) : 0;
    params->virtualTexturePageSize = static_cast<int>(config_.pageSize);
    params->virtualTextureStreamPixels = 0ull;
    params->virtualTextureStreamPageCountX = finalized_ ? static_cast<int>(streamPageCountX_) : 0;
    params->virtualTextureStreamPageCountY = finalized_ ? static_cast<int>(streamPageCountY_) : 0;
    params->virtualTexturePhysicalTextures =
        finalized_ ? reinterpret_cast<unsigned long long>(streamTextureHandlesDevice_.data()) : 0ull;
    params->virtualTexturePhysicalTextureCount = finalized_ ? static_cast<int>(streamPageCountY_) : 0;
    params->virtualTexturePhysicalPagesPerTexture =
        finalized_ ? static_cast<int>(streamPageCountX_) : 0;
    params->virtualTexturePhysicalTextureFormat =
        finalized_ ? static_cast<int>(TextureFormat::RGBA8_UNORM) : 0;
    params->virtualTextureSampleMip = 0;
    params->virtualTextureTextureMeta =
        finalized_ ? reinterpret_cast<unsigned long long>(textureMetaDevice_.data()) : 0ull;
    params->virtualTextureTextureMetaCount = finalized_ ? static_cast<int>(textureMetaHost_.size()) : 0;
    params->virtualTextureMipInfos =
        finalized_ ? reinterpret_cast<unsigned long long>(mipInfosDevice_.data()) : 0ull;
    params->virtualTextureMipInfoCount = finalized_ ? static_cast<int>(mipInfosHost_.size()) : 0;
    params->virtualTextureUdimInfos =
        finalized_ ? reinterpret_cast<unsigned long long>(udimInfosDevice_.data()) : 0ull;
    params->virtualTextureUdimInfoCount = finalized_ ? static_cast<int>(udimInfosHost_.size()) : 0;

}

uint32_t VirtualTextureManager::PackPageTableEntry(uint32_t page,
                                                   uint32_t physicalTextureID,
                                                   uint32_t pageType) const
{
    return ((page & 0xffu) << 0u) | ((physicalTextureID & 0x7fffffu) << 8u) | ((pageType & 0x1u) << 31u);
}

void VirtualTextureManager::UnpackPageTableEntry(uint32_t packed,
                                                 uint32_t *outPage,
                                                 uint32_t *outPhysicalTextureID,
                                                 uint32_t *outPageType) const
{
    if (outPage)
    {
        *outPage = (packed >> 0u) & 0xffu;
    }
    if (outPhysicalTextureID)
    {
        *outPhysicalTextureID = (packed >> 8u) & 0x7fffffu;
    }
    if (outPageType)
    {
        *outPageType = (packed >> 31u) & 0x1u;
    }
}

bool VirtualTextureManager::UpdatePageTableTexel(uint32_t mip,
                                                 uint32_t x,
                                                 uint32_t y,
                                                 uint32_t entry,
                                                 std::string *outError)
{
    if (mip >= pageTableMipOffsets_.size())
    {
        return false;
    }
    if (x >= pageTableMipWidths_[mip] || y >= pageTableMipHeights_[mip])
    {
        return false;
    }
    const uint32_t index = pageTableMipOffsets_[mip] + y * pageTableMipWidths_[mip] + x;
    if (index >= pageTableEntriesHost_.size())
    {
        return false;
    }
    pageTableEntriesHost_[index] = entry;
    if (finalized_ && pageTableEntriesDevice_.data() != nullptr)
    {
        const size_t byteOffset = static_cast<size_t>(index) * sizeof(uint32_t);
        DeviceMemoryView<uint8_t> dst = {pageTableEntriesDevice_.data() + byteOffset, sizeof(uint32_t)};
        device_->CopyBytesToDevice(dst, &entry, sizeof(uint32_t));
    }
    (void)outError;
    return true;
}

bool VirtualTextureManager::AllocateStreamTexture(std::string *outError)
{
    if (!device_)
    {
        if (outError)
        {
            *outError = "VirtualTextureManager: device not initialized";
        }
        return false;
    }
    if (streamTextures_.size() >= maxStreamTextureCount_)
    {
        if (outError)
        {
            *outError = "VirtualTextureManager: stream texture pool exhausted";
        }
        return false;
    }

    const uint32_t physicalTextureID = static_cast<uint32_t>(streamTextures_.size());
    const uint32_t textureWidth = config_.pageSize * streamPageCountX_;
    const uint32_t textureHeight = config_.pageSize;
    std::vector<uint8_t> zeros(streamTextureBytes_, 0u);
    DeviceTextureCreateInfo createInfo = {};
    createInfo.pixels = zeros.data();
    createInfo.pixelBytes = zeros.size();
    createInfo.width = textureWidth;
    createInfo.height = textureHeight;
    createInfo.wrapS = DeviceTextureWrapMode::Clamp;
    createInfo.wrapT = DeviceTextureWrapMode::Clamp;
    createInfo.filter = DeviceTextureFilterMode::Nearest;
    createInfo.format = TextureFormat::RGBA8_UNORM;

    DeviceTexture texture = {};
    if (!device_->CreateTexture(createInfo, &texture, outError))
    {
        return false;
    }

    streamTextures_.push_back(texture);
    if (physicalTextureID < streamTextureHandlesHost_.size())
    {
        streamTextureHandlesHost_[physicalTextureID] = texture.handle;
        if (streamTextureHandlesDevice_.data() != nullptr)
        {
            const size_t byteOffset = static_cast<size_t>(physicalTextureID) * sizeof(unsigned long long);
            DeviceMemoryView<uint8_t> dst = {
                streamTextureHandlesDevice_.data() + byteOffset, sizeof(unsigned long long)};
            device_->CopyBytesToDevice(dst, &streamTextureHandlesHost_[physicalTextureID], sizeof(unsigned long long));
        }
    }
    streamPageCountY_ = static_cast<uint32_t>(streamTextures_.size());

    const uint32_t slotBase = physicalTextureID * streamPageCountX_;
    streamSlots_.resize(static_cast<size_t>(slotBase) + static_cast<size_t>(streamPageCountX_));
    for (uint32_t page = 0u; page < streamPageCountX_; ++page)
    {
        freeSlots_.push_back(slotBase + page);
    }
    streamPageCount_ = static_cast<uint32_t>(streamSlots_.size());
    return true;
}

void VirtualTextureManager::Shutdown()
{
    if (device_)
    {
        for (TextureState &texture : textures_)
        {
            if (texture.tailPixelsDevice.data() != nullptr)
            {
                device_->FreeBytes(texture.tailPixelsDevice);
            }
        }
        for (DeviceTexture &texture : streamTextures_)
        {
            if (texture.valid)
            {
                device_->DestroyTexture(texture);
            }
        }
        if (streamTextureHandlesDevice_.data() != nullptr)
        {
            device_->FreeBytes(streamTextureHandlesDevice_);
        }
        if (textureMetaDevice_.data() != nullptr)
        {
            device_->FreeBytes(textureMetaDevice_);
        }
        if (mipInfosDevice_.data() != nullptr)
        {
            device_->FreeBytes(mipInfosDevice_);
        }
        if (udimInfosDevice_.data() != nullptr)
        {
            device_->FreeBytes(udimInfosDevice_);
        }
        if (pageTableMipHeightsDevice_.data() != nullptr)
        {
            device_->FreeBytes(pageTableMipHeightsDevice_);
        }
        if (pageTableMipWidthsDevice_.data() != nullptr)
        {
            device_->FreeBytes(pageTableMipWidthsDevice_);
        }
        if (pageTableMipOffsetsDevice_.data() != nullptr)
        {
            device_->FreeBytes(pageTableMipOffsetsDevice_);
        }
        if (pageTableEntriesDevice_.data() != nullptr)
        {
            device_->FreeBytes(pageTableEntriesDevice_);
        }
    }

    device_ = nullptr;
    finalized_ = false;
    textures_.clear();
    textureIdToIndex_.clear();
    pageTableMipOffsets_.clear();
    pageTableMipWidths_.clear();
    pageTableMipHeights_.clear();
    pageTableEntriesHost_.clear();
    mipInfosHost_.clear();
    udimInfosHost_.clear();
    textureMetaHost_.clear();
    streamTextures_.clear();
    streamTextureHandlesHost_.clear();
    streamSlots_.clear();
    keyToStreamSlot_.clear();
    lruSlots_.clear();
    freeSlots_.clear();
    streamPageCountX_ = 0u;
    streamPageCountY_ = 0u;
    streamPageCount_ = 0u;
    maxStreamTextureCount_ = 0u;
    streamTextureBytes_ = 0u;
    streamTextureHandlesDevice_ = {};
}

} // namespace texture
} // namespace ybi
