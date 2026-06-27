#pragma once

#include "core/assert.h"
#include <cstdint>
#include <cstdio>
#include <string>

#if defined(_WIN32)

#include <windows.h>

#else
#include <fcntl.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <unistd.h>
#endif

namespace ybi
{

#if defined(_WIN32)
static inline bool OS_UnmapFile(void *ptr)
{
    if (!ptr)
    {
        return true;
    }

    return UnmapViewOfFile(ptr) != 0;
}

static inline bool OS_UnmapFile(void *ptr, uint64_t size)
{
    (void)size;
    return OS_UnmapFile(ptr);
}

static inline bool OS_FlushMappedFile(void *ptr, uint64_t size)
{
    if (!ptr || size == 0)
    {
        return true;
    }

    return FlushViewOfFile(ptr, static_cast<SIZE_T>(size)) != 0;
}

static inline bool OS_ResizeFile(const std::string &filename, uint64_t size)
{
    HANDLE file = CreateFileA(filename.c_str(),
                              GENERIC_READ | GENERIC_WRITE,
                              0,
                              nullptr,
                              OPEN_ALWAYS,
                              FILE_ATTRIBUTE_NORMAL,
                              nullptr);
    if (file == INVALID_HANDLE_VALUE)
    {
        fprintf(stderr, "CreateFileA failed for %s: %lu\n", filename.c_str(), GetLastError());
        return false;
    }

    LARGE_INTEGER newFileSize = {};
    newFileSize.QuadPart = static_cast<LONGLONG>(size);

    const bool success =
        SetFilePointerEx(file, newFileSize, nullptr, FILE_BEGIN) != 0 && SetEndOfFile(file) != 0;
    if (!success)
    {
        fprintf(stderr, "Failed to resize file %s: %lu\n", filename.c_str(), GetLastError());
    }

    CloseHandle(file);
    return success;
}

static inline uint8_t *OS_MapFileWrite(const std::string &filename, uint64_t size)
{
    YBI_ASSERT(size > 0);

    HANDLE file = CreateFileA(filename.c_str(),
                              GENERIC_READ | GENERIC_WRITE,
                              0,
                              nullptr,
                              CREATE_ALWAYS,
                              FILE_ATTRIBUTE_NORMAL,
                              nullptr);
    if (file == INVALID_HANDLE_VALUE)
    {
        fprintf(stderr, "CreateFileA failed for %s: %lu\n", filename.c_str(), GetLastError());
        return nullptr;
    }

    LARGE_INTEGER newFileSize = {};
    newFileSize.QuadPart = static_cast<LONGLONG>(size);
    if (SetFilePointerEx(file, newFileSize, nullptr, FILE_BEGIN) == 0 || SetEndOfFile(file) == 0)
    {
        fprintf(stderr, "Failed to resize file %s: %lu\n", filename.c_str(), GetLastError());
        CloseHandle(file);
        return nullptr;
    }

    HANDLE mapping = CreateFileMappingA(
        file, nullptr, PAGE_READWRITE, newFileSize.HighPart, newFileSize.LowPart, nullptr);
    CloseHandle(file);
    if (!mapping)
    {
        fprintf(
            stderr, "CreateFileMappingA failed for %s: %lu\n", filename.c_str(), GetLastError());
        return nullptr;
    }

    void *ptr = MapViewOfFile(mapping, FILE_MAP_ALL_ACCESS, 0, 0, 0);
    CloseHandle(mapping);
    if (!ptr)
    {
        fprintf(stderr, "MapViewOfFile failed for %s: %lu\n", filename.c_str(), GetLastError());
    }

    return static_cast<uint8_t *>(ptr);
}

static inline uint8_t *OS_MapFileAppend(const std::string &filename, uint64_t size)
{
    YBI_ASSERT(size > 0);

    HANDLE file = CreateFileA(filename.c_str(),
                              GENERIC_READ | GENERIC_WRITE,
                              0,
                              nullptr,
                              OPEN_EXISTING,
                              FILE_ATTRIBUTE_NORMAL,
                              nullptr);
    if (file == INVALID_HANDLE_VALUE)
    {
        fprintf(stderr, "CreateFileA failed for %s: %lu\n", filename.c_str(), GetLastError());
        return nullptr;
    }

    LARGE_INTEGER currentSize = {};
    if (GetFileSizeEx(file, &currentSize) == 0)
    {
        fprintf(stderr, "GetFileSizeEx failed for %s: %lu\n", filename.c_str(), GetLastError());
        CloseHandle(file);
        return nullptr;
    }

    LARGE_INTEGER newFileSize = {};
    newFileSize.QuadPart = currentSize.QuadPart + static_cast<LONGLONG>(size);
    if (SetFilePointerEx(file, newFileSize, nullptr, FILE_BEGIN) == 0 || SetEndOfFile(file) == 0)
    {
        fprintf(stderr, "Failed to resize file %s: %lu\n", filename.c_str(), GetLastError());
        CloseHandle(file);
        return nullptr;
    }

    HANDLE mapping = CreateFileMappingA(
        file, nullptr, PAGE_READWRITE, newFileSize.HighPart, newFileSize.LowPart, nullptr);
    CloseHandle(file);
    if (!mapping)
    {
        fprintf(
            stderr, "CreateFileMappingA failed for %s: %lu\n", filename.c_str(), GetLastError());
        return nullptr;
    }

    void *ptr = MapViewOfFile(mapping, FILE_MAP_ALL_ACCESS, 0, 0, 0);
    CloseHandle(mapping);
    if (!ptr)
    {
        fprintf(stderr, "MapViewOfFile failed for %s: %lu\n", filename.c_str(), GetLastError());
    }

    return static_cast<uint8_t *>(ptr);
}

#else
static inline bool OS_UnmapFile(void *ptr)
{
    if (!ptr)
    {
        return true;
    }

    YBI_ASSERT(false);
    return false;
}

static inline bool OS_UnmapFile(void *ptr, uint64_t size)
{
    if (!ptr)
    {
        return true;
    }

    YBI_ASSERT(size > 0);
    return munmap(ptr, static_cast<size_t>(size)) == 0;
}

static inline bool OS_FlushMappedFile(void *ptr, uint64_t size)
{
    if (!ptr || size == 0)
    {
        return true;
    }

    return msync(ptr, static_cast<size_t>(size), MS_SYNC) == 0;
}

static inline bool OS_ResizeFile(const std::string &filename, uint64_t size)
{
    int file = open(filename.c_str(), O_RDWR | O_CREAT, 0666);
    if (file < 0)
    {
        perror(filename.c_str());
        return false;
    }

    const bool success = ftruncate(file, static_cast<off_t>(size)) == 0;
    if (!success)
    {
        perror(filename.c_str());
    }

    close(file);
    return success;
}

static inline uint8_t *OS_MapFileWrite(const std::string &filename, uint64_t size)
{
    YBI_ASSERT(size > 0);

    int file = open(filename.c_str(), O_RDWR | O_CREAT | O_TRUNC, 0666);
    if (file < 0)
    {
        perror(filename.c_str());
        return nullptr;
    }

    if (ftruncate(file, static_cast<off_t>(size)) != 0)
    {
        perror(filename.c_str());
        close(file);
        return nullptr;
    }

    void *ptr =
        mmap(nullptr, static_cast<size_t>(size), PROT_READ | PROT_WRITE, MAP_SHARED, file, 0);
    close(file);
    if (ptr == MAP_FAILED)
    {
        perror(filename.c_str());
        return nullptr;
    }

    return static_cast<uint8_t *>(ptr);
}

static inline uint8_t *OS_MapFileAppend(const std::string &filename, uint64_t size)
{
    YBI_ASSERT(size > 0);

    int file = open(filename.c_str(), O_RDWR);
    if (file < 0)
    {
        perror(filename.c_str());
        return nullptr;
    }

    struct stat fileStat = {};
    if (fstat(file, &fileStat) != 0)
    {
        perror(filename.c_str());
        close(file);
        return nullptr;
    }

    const uint64_t newFileSize = static_cast<uint64_t>(fileStat.st_size) + size;
    if (ftruncate(file, static_cast<off_t>(newFileSize)) != 0)
    {
        perror(filename.c_str());
        close(file);
        return nullptr;
    }

    void *ptr = mmap(
        nullptr, static_cast<size_t>(newFileSize), PROT_READ | PROT_WRITE, MAP_SHARED, file, 0);
    close(file);
    if (ptr == MAP_FAILED)
    {
        perror(filename.c_str());
        return nullptr;
    }

    return static_cast<uint8_t *>(ptr);
}

#endif
} // namespace ybi
