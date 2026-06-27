cmake_minimum_required(VERSION 3.16)

project(Yeoubi LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17 CACHE STRING "C++ standard for yeoubi and subprojects" FORCE)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/src/cmake")

# Compile Options
option(WITH_CUDA "Enable CUDA backend" OFF)
option(WITH_VULKAN "Enable Vulkan presentation backend" OFF)

set(YBI_COMPILE_DEFINITIONS)
set(YBI_DEPENDENCIES)
set(YBI_SOURCE)
if (WIN32)
    list(APPEND YBI_COMPILE_DEFINITIONS NOMINMAX)
    list(APPEND YBI_COMPILE_DEFINITIONS WIN32_LEAN_AND_MEAN)
endif()

# CUDA
if (WITH_CUDA)
    include(CheckLanguage)
    check_language(CUDA)
    set(YBI_HAS_CUDA OFF)
    if (CMAKE_CUDA_COMPILER)
        set(YBI_HAS_CUDA ON)
    endif()

    if (YBI_HAS_CUDA)
        set(YBI_MIN_OPTIX_VERSION 70300)
        find_package(CUDAToolkit REQUIRED)
        enable_language(CUDA)
        message(STATUS "CUDA Found: ${CMAKE_CUDA_COMPILER} (${CUDA_VERSION})")
        list (APPEND YBI_DEPENDENCIES CUDA::cuda_driver CUDA::cudart)
        list (APPEND YBI_COMPILE_DEFINITIONS WITH_CUDA)

        find_package(OptiX)
        if (OPTIX_FOUND)
            if (OPTIX_VERSION VERSION_LESS ${YBI_MIN_OPTIX_VERSION})
                message(FATAL_ERROR
                    "OptiX ${YBI_MIN_OPTIX_VERSION}+ is required. Found ${OPTIX_VERSION} at ${OPTIX_INCLUDE_DIR}")
            endif()
            list (APPEND YBI_COMPILE_DEFINITIONS WITH_OPTIX OPTIX_VERSION=${OPTIX_VERSION})
            list (APPEND YBI_THIRD_PARTY_INCLUDE_DIRECTORIES "${OPTIX_INCLUDE_DIR}")
            list (APPEND YBI_SOURCE ${CMAKE_CURRENT_SOURCE_DIR}/src/device/cuda_device.cpp)
            # list (APPEND YBI_SOURCE ${CMAKE_CURRENT_SOURCE_DIR}/src/device/bvh/optix/bvh.cpp)
            message (STATUS "OptiX Found: ${OPTIX_INCLUDE_DIR} ${OPTIX_VERSION}")
        else()
            message(FATAL_ERROR "OptiX ${YBI_MIN_OPTIX_VERSION}+ is required if using CUDA.")
        endif()
    endif()
endif()

# Vulkan
if (WITH_VULKAN)
    find_package(Vulkan REQUIRED)
    if (Vulkan_FOUND)
        message(STATUS "Found Vulkan: ${Vulkan_INCLUDE_DIRS}")
    endif()
    list (APPEND YBI_DEPENDENCIES Vulkan::Vulkan)
    list (APPEND YBI_COMPILE_DEFINITIONS WITH_VULKAN)
    list (APPEND YBI_SOURCE ${CMAKE_CURRENT_SOURCE_DIR}/src/display_device/vulkan/vulkan.cpp)
endif()

# GLFW
if (WITH_VULKAN)
    set(GLFW_BUILD_DOCS OFF CACHE BOOL "" FORCE)
    set(GLFW_BUILD_TESTS OFF CACHE BOOL "" FORCE)
    set(GLFW_BUILD_EXAMPLES OFF CACHE BOOL "" FORCE)
    add_subdirectory(third_party/glfw)
    list(APPEND YBI_DEPENDENCIES glfw)
endif()

# USD
set(USD_ROOT "" CACHE PATH "Path to the OpenUSD installation directory")
if (USD_ROOT)
    find_package(OpenUSD)
    if (OPENUSD_FOUND)
        set (WITH_USD ON)
        set (OPENUSD_INCLUDE_DIRS ${OPENUSD_INCLUDE_DIRS})
        set (OPENUSD_LIBRARIES ${OPENUSD_LIBRARIES})
        set (OPENUSD_LIBRARY_DIR ${OPENUSD_LIBRARY_DIR})
    endif()
    list (APPEND YBI_DEPENDENCIES ${OPENUSD_LIBRARIES})
endif()

# Virtual textures exe
add_executable(yeoubi_virtual_textures
               ${CMAKE_CURRENT_SOURCE_DIR}/src/core/allocation.cpp
               ${CMAKE_CURRENT_SOURCE_DIR}/src/scene/scene.cpp
               ${CMAKE_CURRENT_SOURCE_DIR}/src/cmd/virtual_textures.cpp
               ${CMAKE_CURRENT_SOURCE_DIR}/src/core/string_builder.cpp
)

target_include_directories(yeoubi_virtual_textures PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/src
    ${OPENUSD_INCLUDE_DIRS}
)

target_include_directories(yeoubi_virtual_textures SYSTEM PRIVATE
    "${CMAKE_CURRENT_SOURCE_DIR}/third_party"
)

target_link_libraries(yeoubi_virtual_textures PRIVATE
    ${OPENUSD_LIBRARIES}
)

target_compile_definitions(yeoubi_virtual_textures PRIVATE
    ${YBI_COMPILE_DEFINITIONS}
    $<$<CONFIG:Debug,RelWithDebInfo>:YBI_DEBUG>
)

# Renderer exe
add_executable(yeoubi
               ${CMAKE_CURRENT_SOURCE_DIR}/src/core/allocation.cpp
               ${CMAKE_CURRENT_SOURCE_DIR}/src/io/usd.cpp
               ${CMAKE_CURRENT_SOURCE_DIR}/src/scene/scene.cpp
               ${CMAKE_CURRENT_SOURCE_DIR}/src/device/device.cpp
               ${CMAKE_CURRENT_SOURCE_DIR}/src/display_device/display_device.cpp
               ${CMAKE_CURRENT_SOURCE_DIR}/src/input/camera_controller.cpp
               ${CMAKE_CURRENT_SOURCE_DIR}/src/main.cpp
               ${YBI_SOURCE}
)

set (YBI_KERNEL_HEADERS
    ${CMAKE_CURRENT_SOURCE_DIR}/src/core/vec3f.h
)

list(APPEND YBI_THIRD_PARTY_INCLUDE_DIRECTORIES
    "${CMAKE_CURRENT_SOURCE_DIR}/third_party"
)


if (YBI_HAS_CUDA AND OPTIX_FOUND)
    set(YBI_OPTIX_PTX_PATH "${CMAKE_CURRENT_BINARY_DIR}/kernel.ptx")
    add_custom_command(
        OUTPUT "${YBI_OPTIX_PTX_PATH}"
        COMMAND "${CMAKE_CUDA_COMPILER}"
        --ptx
        "--gpu-architecture=compute_75"
        "--std=c++17"
        # "--expt-relaxed-constexpr"
        "--use_fast_math"
        "-I${CMAKE_CURRENT_SOURCE_DIR}/src"
        "-I${OPTIX_INCLUDE_DIR}"
        # "-I${YBI_CUDA_INCLUDE_DIR}"
        # "-I${OPENPBR_BSDF_INCLUDE_DIR}"
        "-o" "${YBI_OPTIX_PTX_PATH}"
        "${CMAKE_CURRENT_SOURCE_DIR}/src/render/optix/kernel.cu"
        DEPENDS
            "${CMAKE_CURRENT_SOURCE_DIR}/src/render/optix/kernel.cu"
            ${YBI_KERNEL_HEADERS}
        VERBATIM
    )
    add_custom_target(yeoubi_optix_kernel DEPENDS ${YBI_OPTIX_PTX_PATH})
    add_dependencies(yeoubi yeoubi_optix_kernel)

    list (APPEND YBI_COMPILE_DEFINITIONS YBI_OPTIX_PTX_PATH=\"${YBI_OPTIX_PTX_PATH}\")
endif()

target_include_directories(yeoubi PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/src
    ${OPENUSD_INCLUDE_DIRS}
)

target_include_directories(yeoubi SYSTEM PRIVATE
    ${YBI_THIRD_PARTY_INCLUDE_DIRECTORIES}
)

target_link_libraries(yeoubi PRIVATE
    ${YBI_DEPENDENCIES}
)

target_compile_definitions(yeoubi PRIVATE
    ${YBI_COMPILE_DEFINITIONS}
    $<$<CONFIG:Debug,RelWithDebInfo>:YBI_DEBUG>
)
