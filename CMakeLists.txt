cmake_minimum_required(VERSION 3.16)

project(Yeoubi LANGUAGES C CXX)

set(CMAKE_CXX_STANDARD 17 CACHE STRING "C++ standard for yeoubi and subprojects" FORCE)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/src/cmake")

# SIMD feature detection
include(CheckCXXSourceCompiles)
check_cxx_source_compiles(
    "#include <xmmintrin.h>
     int main()
     {
         const __m128 result = _mm_rsqrt_ss(_mm_set_ss(1.0f));
         return _mm_cvtss_f32(result) == 0.0f;
     }"
    YBI_HAS_SSE)

# Compile Options
option(WITH_CUDA "Enable CUDA backend" OFF)
option(WITH_VULKAN "Enable Vulkan presentation backend" OFF)
option(WITH_PATH_GUIDING "Enable OpenPGL path guiding" OFF)

set(YBI_COMPILE_DEFINITIONS)
set(YBI_DEPENDENCIES)
set(YBI_SOURCE)
set(YBI_USE_GPU)
if (YBI_HAS_SSE)
    list(APPEND YBI_COMPILE_DEFINITIONS YBI_HAS_SSE)
endif()
if (WIN32)
    list(APPEND YBI_COMPILE_DEFINITIONS NOMINMAX)
    list(APPEND YBI_COMPILE_DEFINITIONS WIN32_LEAN_AND_MEAN)
endif()

# CUDA
if (WITH_CUDA)
    set(YBI_USE_GPU ON)
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
            list (APPEND YBI_SOURCE ${CMAKE_CURRENT_SOURCE_DIR}/src/device/cuda/clusters.cu)
            list (APPEND YBI_SOURCE
                ${CMAKE_CURRENT_SOURCE_DIR}/src/device/cuda/validate_device_cluster.cu)
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

# meshoptimizer
set(MESHOPT_BUILD_DEMO OFF CACHE BOOL "" FORCE)
set(MESHOPT_BUILD_GLTFPACK OFF CACHE BOOL "" FORCE)
set(MESHOPT_BUILD_SHARED_LIBS OFF CACHE BOOL "" FORCE)
set(MESHOPT_INSTALL OFF CACHE BOOL "" FORCE)
add_subdirectory(third_party/meshoptimizer EXCLUDE_FROM_ALL)
list(APPEND YBI_DEPENDENCIES meshoptimizer)

list(APPEND YBI_THIRD_PARTY_INCLUDE_DIRECTORIES "${CMAKE_CURRENT_SOURCE_DIR}/third_party/meshoptimizer/src")

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

# OpenPGL
if (WITH_PATH_GUIDING)
    find_package(openpgl CONFIG REQUIRED)
    if (openpgl_FOUND)
        list (APPEND YBI_DEPENDENCIES openpgl::openpgl)
        list (APPEND YBI_COMPILE_DEFINITIONS WITH_PATH_GUIDING)
    endif()
endif()

# Virtual textures exe
set(YBI_MINIZ_SOURCE "${CMAKE_CURRENT_SOURCE_DIR}/third_party/tinyexr/deps/miniz/miniz.c")

# LibTIFF
set(tiff-tools OFF CACHE BOOL "" FORCE)
set(tiff-tests OFF CACHE BOOL "" FORCE)
set(tiff-contrib OFF CACHE BOOL "" FORCE)
set(tiff-docs OFF CACHE BOOL "" FORCE)
set(tiff-install OFF CACHE BOOL "" FORCE)
set(tiff-cxx OFF CACHE BOOL "" FORCE)
add_subdirectory(third_party/libtiff EXCLUDE_FROM_ALL)

add_executable(yeoubi_virtual_textures
               ${CMAKE_CURRENT_SOURCE_DIR}/src/core/allocation.cpp
               ${CMAKE_CURRENT_SOURCE_DIR}/src/io/usd.cpp
               ${CMAKE_CURRENT_SOURCE_DIR}/src/io/usd_material.cpp
               ${CMAKE_CURRENT_SOURCE_DIR}/src/scene/shader_compilation.cpp
               ${CMAKE_CURRENT_SOURCE_DIR}/src/scene/scene.cpp
               ${CMAKE_CURRENT_SOURCE_DIR}/src/cmd/virtual_textures.cpp
               ${CMAKE_CURRENT_SOURCE_DIR}/src/core/string_builder.cpp
               ${YBI_MINIZ_SOURCE}
)

target_include_directories(yeoubi_virtual_textures PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/src
    ${OPENUSD_INCLUDE_DIRS}
)

target_include_directories(yeoubi_virtual_textures SYSTEM PRIVATE
    "${CMAKE_CURRENT_SOURCE_DIR}/third_party"
    "${CMAKE_CURRENT_SOURCE_DIR}/third_party/tinyexr/deps/miniz"
)

target_link_libraries(yeoubi_virtual_textures PRIVATE
    ${OPENUSD_LIBRARIES}
    TIFF::tiff
)

target_compile_definitions(yeoubi_virtual_textures PRIVATE
    ${YBI_COMPILE_DEFINITIONS}
    $<$<CONFIG:Debug,RelWithDebInfo>:YBI_DEBUG>
)

# Generate pmj tables exe
add_executable(yeoubi_pmj_tables
                ${CMAKE_CURRENT_SOURCE_DIR}/src/cmd/generate_pmj02_tables.cpp
                ${CMAKE_CURRENT_SOURCE_DIR}/src/render/pmj02.cpp
)

target_include_directories(yeoubi_pmj_tables PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/src
)

target_compile_definitions(yeoubi_pmj_tables PRIVATE
    ${YBI_COMPILE_DEFINITIONS}
    $<$<CONFIG:Debug,RelWithDebInfo>:YBI_DEBUG>
)

# Renderer exe
add_executable(yeoubi
               ${CMAKE_CURRENT_SOURCE_DIR}/src/core/allocation.cpp
               ${CMAKE_CURRENT_SOURCE_DIR}/src/io/usd.cpp
               ${CMAKE_CURRENT_SOURCE_DIR}/src/io/usd_material.cpp
               ${CMAKE_CURRENT_SOURCE_DIR}/src/scene/shader_compilation.cpp
               ${CMAKE_CURRENT_SOURCE_DIR}/src/scene/scene.cpp
               ${CMAKE_CURRENT_SOURCE_DIR}/src/scene/render_system.cpp
               ${CMAKE_CURRENT_SOURCE_DIR}/src/device/device.cpp
               ${CMAKE_CURRENT_SOURCE_DIR}/src/display_device/display_device.cpp
               ${CMAKE_CURRENT_SOURCE_DIR}/src/input/camera_controller.cpp
               ${CMAKE_CURRENT_SOURCE_DIR}/src/render/virtual_texture_manager.cpp
               ${CMAKE_CURRENT_SOURCE_DIR}/src/scene/clusterizer.cpp
               ${CMAKE_CURRENT_SOURCE_DIR}/src/scene/lut.cpp
               ${CMAKE_CURRENT_SOURCE_DIR}/src/main.cpp
               ${YBI_SOURCE}
)

list(APPEND YBI_THIRD_PARTY_INCLUDE_DIRECTORIES
    "${CMAKE_CURRENT_SOURCE_DIR}/third_party"
)


if (YBI_HAS_CUDA AND OPTIX_FOUND)
    get_filename_component(YBI_CUDA_HOST_COMPILER_DIR "${CMAKE_CXX_COMPILER}" DIRECTORY)

    set(YBI_OPTIX_PTX_PATH "${CMAKE_CURRENT_BINARY_DIR}/kernel.ptx")
    set(YBI_OPTIX_PTX_DEPFILE "${CMAKE_CURRENT_BINARY_DIR}/kernel.ptx.d")
    add_custom_command(
        OUTPUT "${YBI_OPTIX_PTX_PATH}"
        COMMAND "${CMAKE_CUDA_COMPILER}"
        --ptx
        "--gpu-architecture=compute_75"
        "--std=c++17"
        "--compiler-bindir=${YBI_CUDA_HOST_COMPILER_DIR}"
        # "--expt-relaxed-constexpr"
        "--use_fast_math"
        "-DOPTIX_VERSION=${OPTIX_VERSION}"
        "-I${CMAKE_CURRENT_SOURCE_DIR}/src"
        "-I${OPTIX_INCLUDE_DIR}"
        # "-I${YBI_CUDA_INCLUDE_DIR}"
        # "-I${OPENPBR_BSDF_INCLUDE_DIR}"
        -MD
        -MF "${YBI_OPTIX_PTX_DEPFILE}"
        "-o" "${YBI_OPTIX_PTX_PATH}"
        "${CMAKE_CURRENT_SOURCE_DIR}/src/render/optix/kernel.cu"
        DEPENDS
            "${CMAKE_CURRENT_SOURCE_DIR}/src/render/optix/kernel.cu"
        DEPFILE "${YBI_OPTIX_PTX_DEPFILE}"
        VERBATIM
    )
    add_custom_target(yeoubi_optix_kernel DEPENDS ${YBI_OPTIX_PTX_PATH})
    add_dependencies(yeoubi yeoubi_optix_kernel)

    set(YBI_VIRTUAL_TEXTURE_KERNEL_PATH "${CMAKE_CURRENT_BINARY_DIR}/virtual_textures.ptx")
    set(YBI_VIRTUAL_TEXTURE_KERNEL_DEPFILE
        "${CMAKE_CURRENT_BINARY_DIR}/virtual_textures.ptx.d")
    add_custom_command(
        OUTPUT "${YBI_VIRTUAL_TEXTURE_KERNEL_PATH}"
        COMMAND "${CMAKE_CUDA_COMPILER}"
        --ptx
        "--gpu-architecture=sm_75"
        "--std=c++17"
        "--compiler-bindir=${YBI_CUDA_HOST_COMPILER_DIR}"
        "--use_fast_math"
        "-I${CMAKE_CURRENT_SOURCE_DIR}/src"
        "-Ofc=max"
        -MD
        -MF "${YBI_VIRTUAL_TEXTURE_KERNEL_DEPFILE}"
        "-o" "${YBI_VIRTUAL_TEXTURE_KERNEL_PATH}"
        "${CMAKE_CURRENT_SOURCE_DIR}/src/render/cuda/virtual_textures.cu"
        DEPENDS
            "${CMAKE_CURRENT_SOURCE_DIR}/src/render/cuda/virtual_textures.cu"
        DEPFILE "${YBI_VIRTUAL_TEXTURE_KERNEL_DEPFILE}"
        VERBATIM
    )
    add_custom_target(yeoubi_virtual_texture_kernel DEPENDS ${YBI_VIRTUAL_TEXTURE_KERNEL_PATH})
    add_dependencies(yeoubi yeoubi_virtual_texture_kernel)

    list (APPEND YBI_COMPILE_DEFINITIONS YBI_OPTIX_PTX_PATH=\"${YBI_OPTIX_PTX_PATH}\")
    list (APPEND YBI_COMPILE_DEFINITIONS YBI_VIRTUAL_TEXTURE_KERNEL_PATH=\"${YBI_VIRTUAL_TEXTURE_KERNEL_PATH}\")
endif()

target_include_directories(yeoubi PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/src
    ${OPENUSD_INCLUDE_DIRS}
    ${YBI_INCLUDE_DIRS}
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
