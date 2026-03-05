#include "texture/virtual_texture/manager.h"

#include <algorithm>
#include <cmath>

namespace ybi
{
namespace texture
{

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
    uint64_t desiredPages = config_.cacheBytes / std::max<size_t>(1u, pageBytes);
    if (desiredPages == 0u)
    {
        desiredPages = 1u;
    }
    uint32_t pageCountX = static_cast<uint32_t>(std::ceil(std::sqrt(double(desiredPages))));
    pageCountX = std::max(1u, std::min(pageCountX, 4095u));
    uint32_t pageCountY = static_cast<uint32_t>((desiredPages + pageCountX - 1u) / pageCountX);
    pageCountY = std::max(1u, std::min(pageCountY, 4095u));
    streamPageCountX_ = pageCountX;
    streamPageCountY_ = pageCountY;
    streamPageCount_ = streamPageCountX_ * streamPageCountY_;

    streamPixelsHost_.assign(static_cast<size_t>(streamPageCount_) * pageBytes, 0u);
    streamPixelsDevice_ = device_->AllocBytes(streamPixelsHost_.size());
    device_->CopyBytesToDevice(streamPixelsDevice_, streamPixelsHost_.data(), streamPixelsHost_.size());

    streamSlots_.assign(streamPageCount_, StreamSlotState{});
    keyToStreamSlot_.clear();
    lruSlots_.clear();
    freeSlots_.clear();
    freeSlots_.reserve(streamPageCount_);
    for (uint32_t slot = 0u; slot < streamPageCount_; ++slot)
    {
        freeSlots_.push_back(slot);
    }
    (void)outError;
    return true;
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
    params->virtualTextureStreamPixels =
        finalized_ ? reinterpret_cast<unsigned long long>(streamPixelsDevice_.data()) : 0ull;
    params->virtualTextureStreamPageCountX = finalized_ ? static_cast<int>(streamPageCountX_) : 0;
    params->virtualTextureStreamPageCountY = finalized_ ? static_cast<int>(streamPageCountY_) : 0;
    params->virtualTextureSampleMip = 1;
    params->virtualTextureTextureMeta =
        finalized_ ? reinterpret_cast<unsigned long long>(textureMetaDevice_.data()) : 0ull;
    params->virtualTextureTextureMetaCount = finalized_ ? static_cast<int>(textureMetaHost_.size()) : 0;
    params->virtualTextureMipInfos =
        finalized_ ? reinterpret_cast<unsigned long long>(mipInfosDevice_.data()) : 0ull;
    params->virtualTextureMipInfoCount = finalized_ ? static_cast<int>(mipInfosHost_.size()) : 0;

    params->virtualTextureTileEntries = 0ull;
    params->virtualTextureTilePixels = 0ull;
    params->virtualTextureTileEntryCapacity = 0;
}

uint32_t VirtualTextureManager::PackPageTableEntry(uint32_t pageX,
                                                   uint32_t pageY,
                                                   uint32_t pageType,
                                                   uint32_t flags) const
{
    return ((pageX & 0xfffu) << 0u) | ((pageY & 0xfffu) << 12u) | ((pageType & 0xfu) << 24u) |
           ((flags & 0xfu) << 28u);
}

void VirtualTextureManager::UnpackPageTableEntry(uint32_t packed,
                                                 uint32_t *outPageX,
                                                 uint32_t *outPageY,
                                                 uint32_t *outPageType,
                                                 uint32_t *outFlags) const
{
    if (outPageX)
    {
        *outPageX = (packed >> 0u) & 0xfffu;
    }
    if (outPageY)
    {
        *outPageY = (packed >> 12u) & 0xfffu;
    }
    if (outPageType)
    {
        *outPageType = (packed >> 24u) & 0xfu;
    }
    if (outFlags)
    {
        *outFlags = (packed >> 28u) & 0xfu;
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
        if (streamPixelsDevice_.data() != nullptr)
        {
            device_->FreeBytes(streamPixelsDevice_);
        }
        if (textureMetaDevice_.data() != nullptr)
        {
            device_->FreeBytes(textureMetaDevice_);
        }
        if (mipInfosDevice_.data() != nullptr)
        {
            device_->FreeBytes(mipInfosDevice_);
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
    textureMetaHost_.clear();
    streamPixelsHost_.clear();
    streamSlots_.clear();
    keyToStreamSlot_.clear();
    lruSlots_.clear();
    freeSlots_.clear();
    streamPageCountX_ = 0u;
    streamPageCountY_ = 0u;
    streamPageCount_ = 0u;
}

} // namespace texture
} // namespace ybi
