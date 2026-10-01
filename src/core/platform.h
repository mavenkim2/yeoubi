#pragma once

#include <cstddef>
#include <filesystem>
#include <string>
#include <vector>

struct GLFWwindow;

#if defined(_WIN32)
#define GLFW_EXPOSE_NATIVE_WIN32
#include <GLFW/glfw3.h>
#include <GLFW/glfw3native.h>
#include <commdlg.h>

namespace ybi
{

inline std::wstring FileDialogFilter(const std::vector<std::filesystem::path> &filePatterns)
{
    std::wstring filter;
    if (!filePatterns.empty())
    {
        filter += L"Supported files";
        filter.push_back(L'\0');
        for (std::size_t i = 0; i < filePatterns.size(); ++i)
        {
            if (i != 0)
                filter.push_back(L';');
            filter += filePatterns[i].wstring();
        }
        filter.push_back(L'\0');
    }
    filter += L"All files";
    filter.push_back(L'\0');
    filter += L"*.*";
    filter.push_back(L'\0');
    filter.push_back(L'\0');
    return filter;
}

inline std::string OpenFileDialog(GLFWwindow *window,
                                  const std::vector<std::filesystem::path> &filePatterns,
                                  std::string *error)
{
    const std::wstring filter = FileDialogFilter(filePatterns);
    wchar_t selected[32768] = {};
    OPENFILENAMEW dialog = {};
    dialog.lStructSize = sizeof(dialog);
    dialog.hwndOwner = glfwGetWin32Window(window);
    dialog.lpstrFilter = filter.c_str();
    dialog.lpstrFile = selected;
    dialog.nMaxFile = static_cast<DWORD>(sizeof(selected) / sizeof(selected[0]));
    dialog.Flags = OFN_FILEMUSTEXIST | OFN_PATHMUSTEXIST | OFN_NOCHANGEDIR;
    if (!GetOpenFileNameW(&dialog))
    {
        if (CommDlgExtendedError() != 0)
            *error = "The Windows file dialog could not open.";
        return {};
    }

    const int length = WideCharToMultiByte(CP_UTF8, 0, selected, -1, nullptr, 0, nullptr, nullptr);
    if (length <= 0)
    {
        *error = "The selected path could not be converted to UTF-8.";
        return {};
    }
    std::string path(static_cast<std::size_t>(length), '\0');
    if (WideCharToMultiByte(CP_UTF8, 0, selected, -1, path.data(), length, nullptr, nullptr) !=
        length)
    {
        *error = "The selected path could not be converted to UTF-8.";
        return {};
    }
    path.pop_back();
    return path;
}

inline std::string SaveFileDialog(GLFWwindow *window,
                                  const std::vector<std::filesystem::path> &filePatterns,
                                  std::string *error)
{
    const std::wstring filter = FileDialogFilter(filePatterns);
    std::wstring defaultExtension;
    if (!filePatterns.empty())
    {
        defaultExtension = filePatterns.front().extension().wstring();
        if (!defaultExtension.empty() && defaultExtension.front() == L'.')
            defaultExtension.erase(0, 1);
    }

    wchar_t selected[32768] = {};
    OPENFILENAMEW dialog = {};
    dialog.lStructSize = sizeof(dialog);
    dialog.hwndOwner = glfwGetWin32Window(window);
    dialog.lpstrFilter = filter.c_str();
    dialog.lpstrFile = selected;
    dialog.nMaxFile = static_cast<DWORD>(sizeof(selected) / sizeof(selected[0]));
    if (!defaultExtension.empty())
        dialog.lpstrDefExt = defaultExtension.c_str();
    dialog.Flags = OFN_OVERWRITEPROMPT | OFN_PATHMUSTEXIST | OFN_NOCHANGEDIR;
    if (!GetSaveFileNameW(&dialog))
    {
        if (CommDlgExtendedError() != 0)
            *error = "The Windows save dialog could not open.";
        return {};
    }

    const int length = WideCharToMultiByte(CP_UTF8, 0, selected, -1, nullptr, 0, nullptr, nullptr);
    if (length <= 0)
    {
        *error = "The selected path could not be converted to UTF-8.";
        return {};
    }
    std::string path(static_cast<std::size_t>(length), '\0');
    if (WideCharToMultiByte(CP_UTF8, 0, selected, -1, path.data(), length, nullptr, nullptr) !=
        length)
    {
        *error = "The selected path could not be converted to UTF-8.";
        return {};
    }
    path.pop_back();
    return path;
}

} // namespace ybi
#endif
