#include "core/string_builder.h"
#include "core/mmap.h"

#include <algorithm>
#include <cstring>
#include <utility>

namespace ybi
{

static constexpr uint64_t kFileBuilderChunkSize = 512ull * 1024ull * 1024ull;

FileBuilder::FileBuilder(std::string filename)
{
    Open(std::move(filename));
}

FileBuilder::~FileBuilder()
{
    Close();
}

void FileBuilder::Open(std::string newFilename)
{
    Close();

    filename = std::move(newFilename);
    totalSize = 0;
    currentLimit = 0;
    ptr = nullptr;
}

void FileBuilder::Close()
{
    if (ptr)
    {
        OS_FlushMappedFile(ptr, totalSize);
        OS_UnmapFile(ptr, currentLimit);
        ptr = nullptr;
        currentLimit = 0;
    }

    if (!filename.empty())
    {
        OS_ResizeFile(filename, totalSize);
    }
}

bool FileBuilder::Expand(uint64_t size)
{
    if (size == 0)
    {
        return true;
    }

    YBI_ASSERT(!filename.empty());

    if (ptr && totalSize + size <= currentLimit)
    {
        return true;
    }

    uint64_t newLimit = currentLimit == 0 ? kFileBuilderChunkSize : currentLimit;
    const uint64_t requiredSize = totalSize + size;
    while (newLimit < requiredSize)
    {
        newLimit += std::max(kFileBuilderChunkSize, size);
    }

    if (ptr)
    {
        OS_FlushMappedFile(ptr, totalSize);
        OS_UnmapFile(ptr, currentLimit);
        OS_ResizeFile(filename, totalSize);
    }

    ptr = currentLimit == 0 ? OS_MapFileWrite(filename, newLimit)
                            : OS_MapFileAppend(filename, newLimit - totalSize);
    if (!ptr)
    {
        currentLimit = 0;
        return false;
    }

    currentLimit = newLimit;
    return true;
}

uint64_t FileBuilder::Allocate(uint64_t size)
{
    const uint64_t offset = totalSize;
    if (!Expand(size))
    {
        return InvalidOffset;
    }

    totalSize += size;
    return offset;
}

uint8_t *FileBuilder::GetPtr(uint64_t offset)
{
    YBI_ASSERT(ptr);
    YBI_ASSERT(offset <= totalSize);

    if (!ptr || offset > totalSize)
    {
        return nullptr;
    }

    return ptr + offset;
}

const uint8_t *FileBuilder::GetPtr(uint64_t offset) const
{
    YBI_ASSERT(ptr);
    YBI_ASSERT(offset <= totalSize);

    if (!ptr || offset > totalSize)
    {
        return nullptr;
    }

    return ptr + offset;
}

uint64_t FileBuilder::AppendBytes(const void *data, size_t size)
{
    YBI_ASSERT(data || size == 0);
    if (size == 0)
    {
        return totalSize;
    }

    const uint64_t offset = Allocate(size);
    if (offset == InvalidOffset)
    {
        return InvalidOffset;
    }

    std::memcpy(GetPtr(offset), data, size);
    return offset;
}

} // namespace ybi
