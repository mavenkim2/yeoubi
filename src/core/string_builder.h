#pragma once

#include "core/arena.h"
#include <cstddef>
#include <cstdint>
#include <string>

namespace ybi
{

class FileBuilder
{
    std::string filename;
    uint8_t *ptr = nullptr;
    uint64_t totalSize = 0;
    uint64_t currentLimit = 0;

    Arena arena;

    bool Expand(uint64_t size);

public:
    static constexpr uint64_t InvalidOffset = ~uint64_t(0);

    uint64_t Allocate(uint64_t size);
    uint8_t *GetPtr(uint64_t offset);
    const uint8_t *GetPtr(uint64_t offset) const;
    uint64_t AppendBytes(const void *data, size_t size);

    FileBuilder() = default;
    explicit FileBuilder(std::string filename);
    ~FileBuilder();

    void Open(std::string filename);
    void Close();

    template <typename T>
    uint64_t Append(const T *data, size_t count)
    {
        return AppendBytes(data, sizeof(T) * count);
    }
};

} // namespace ybi
