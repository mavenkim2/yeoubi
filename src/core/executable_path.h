#pragma once

#include <array>
#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <string>

#if defined(_WIN32)
#include <windows.h>
#elif defined(__linux__)
#include <unistd.h>
#elif defined(__APPLE__)
#include <mach-o/dyld.h>
#include <system_error>
#else
#error Unsupported platform for ExecutableDirectory
#endif

namespace ybi
{

// Return an empty path if the executable path cannot be read or exceeds this limit.
inline std::filesystem::path ExecutableDirectory()
{
    constexpr std::size_t MaxExecutablePathLength = 32768;

#if defined(_WIN32)
    std::array<wchar_t, MaxExecutablePathLength> path{};
    const DWORD length = GetModuleFileNameW(nullptr, path.data(), DWORD(path.size()));
    if (length == 0 || length >= path.size())
        return {};
    return std::filesystem::path(path.data()).parent_path();
#elif defined(__linux__)
    std::array<char, MaxExecutablePathLength> path{};
    const ssize_t length = readlink("/proc/self/exe", path.data(), path.size());
    if (length < 0 || std::size_t(length) >= path.size())
        return {};
    return std::filesystem::path(std::string(path.data(), std::size_t(length))).parent_path();
#else
    std::array<char, MaxExecutablePathLength> path{};
    std::uint32_t size = std::uint32_t(path.size());
    if (_NSGetExecutablePath(path.data(), &size) != 0)
        return {};
    std::error_code error;
    const std::filesystem::path executable = std::filesystem::canonical(path.data(), error);
    return error ? std::filesystem::path{} : executable.parent_path();
#endif
}

} // namespace ybi
