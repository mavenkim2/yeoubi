cmake_minimum_required(VERSION 3.16)

project(Yeoubi LANGUAGES C CXX)

set(CMAKE_CXX_STANDARD 17 CACHE STRING "C++ standard for yeoubi and subprojects" FORCE)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/src/cmake")

# SIMD feature detection
include(CheckCXXSourceCompiles)
include(CMakePushCheckState)
cmake_push_check_state(RESET)
if (CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang" AND NOT MSVC)
    set(CMAKE_REQUIRED_FLAGS "-msse4.1")
endif()
check_cxx_source_compiles(
    "#include <smmintrin.h>
     int main()
     {
         const __m128 result = _mm_blendv_ps(_mm_setzero_ps(), _mm_set1_ps(1.0f),
                                           _mm_castsi128_ps(_mm_set1_epi32(-1)));
         return _mm_cvtss_f32(result) == 0.0f;
     }"
    YBI_HAS_SSE41)
cmake_pop_check_state()

# Compile Options
option(WITH_CUDA "Enable CUDA backend" OFF)
option(WITH_VULKAN "Enable Vulkan presentation backend" OFF)
option(WITH_PATH_GUIDING "Enable OpenPGL path guiding" OFF)
option(YBI_USE_PATH_GUIDING "Enable Yeoubi path guiding" OFF)
option(YBI_MATERIAL_DEBUG "Print MaterialX BSDF sampling decisions" OFF)

set(YBI_COMPILE_DEFINITIONS)
set(YBI_COMPILE_OPTIONS)
set(YBI_DEPENDENCIES)
set(YBI_SOURCE)
set(YBI_USE_GPU)
if (YBI_HAS_SSE41)
    list(APPEND YBI_COMPILE_DEFINITIONS $<$<COMPILE_LANGUAGE:CXX>:SIMD_4>)
    if (CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang" AND NOT MSVC)
        list(APPEND YBI_COMPILE_OPTIONS $<$<COMPILE_LANGUAGE:CXX>:-msse4.1>)
    endif()
endif()
if (YBI_MATERIAL_DEBUG)
    list(APPEND YBI_COMPILE_DEFINITIONS YBI_MATERIAL_DEBUG)
endif()
if (WIN32)
    list(APPEND YBI_COMPILE_DEFINITIONS NOMINMAX)
    list(APPEND YBI_COMPILE_DEFINITIONS WIN32_LEAN_AND_MEAN)
    list(APPEND YBI_DEPENDENCIES comdlg32)
endif()

set (YBI_MATERIALX_STUBS
     ${CMAKE_CURRENT_SOURCE_DIR}/src/kernel/bsdf/materialx/stubs/artistic_ior_extinction_impl.h
     ${CMAKE_CURRENT_SOURCE_DIR}/src/kernel/bsdf/materialx/stubs/artistic_ior_ior_impl.h
     ${CMAKE_CURRENT_SOURCE_DIR}/src/kernel/bsdf/materialx/stubs/closure_impl.h
     ${CMAKE_CURRENT_SOURCE_DIR}/src/kernel/bsdf/materialx/stubs/conductor_impl.h
     ${CMAKE_CURRENT_SOURCE_DIR}/src/kernel/bsdf/materialx/stubs/dielectric_impl.h
     ${CMAKE_CURRENT_SOURCE_DIR}/src/kernel/bsdf/materialx/stubs/image_impl.h
     ${CMAKE_CURRENT_SOURCE_DIR}/src/kernel/bsdf/materialx/stubs/layer_impl.h
     ${CMAKE_CURRENT_SOURCE_DIR}/src/kernel/bsdf/materialx/stubs/mix_impl.h
     ${CMAKE_CURRENT_SOURCE_DIR}/src/kernel/bsdf/materialx/stubs/normalmap_impl.h
     ${CMAKE_CURRENT_SOURCE_DIR}/src/kernel/bsdf/materialx/stubs/oren_nayar_impl.h
     ${CMAKE_CURRENT_SOURCE_DIR}/src/kernel/bsdf/materialx/stubs/remap_impl.h
     ${CMAKE_CURRENT_SOURCE_DIR}/src/kernel/bsdf/materialx/stubs/sheen_impl.h
     ${CMAKE_CURRENT_SOURCE_DIR}/src/kernel/bsdf/materialx/stubs/tiledimage_impl.h
     ${CMAKE_CURRENT_SOURCE_DIR}/src/kernel/bsdf/materialx/stubs/value_impl.h)

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
        list (APPEND YBI_DEPENDENCIES CUDA::cuda_driver CUDA::cudart CUDA::nvrtc)
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
            list (APPEND YBI_SOURCE ${CMAKE_CURRENT_SOURCE_DIR}/src/kernel/backend/cuda/clusters.cu)
            list (APPEND YBI_SOURCE
                ${CMAKE_CURRENT_SOURCE_DIR}/src/kernel/backend/cuda/validate_device_cluster.cu)
            # list (APPEND YBI_SOURCE ${CMAKE_CURRENT_SOURCE_DIR}/src/device/bvh/optix/bvh.cpp)
            message (STATUS "OptiX Found: ${OPTIX_INCLUDE_DIR} ${OPTIX_VERSION}")
        else()
            message(FATAL_ERROR "OptiX ${YBI_MIN_OPTIX_VERSION}+ is required if using CUDA.")
        endif()
    endif()
endif()

if (YBI_USE_PATH_GUIDING)
    if (NOT YBI_HAS_SSE41)
        message(FATAL_ERROR "YBI_USE_PATH_GUIDING requires SSE4.1 for the CPU kd-tree backend")
    endif()
    set(YBI_KD_TREE_KERNEL_PATH "")
    if (YBI_HAS_CUDA AND OPTIX_FOUND)
        set(YBI_KD_TREE_KERNEL_PATH "${CMAKE_CURRENT_BINARY_DIR}/kd_tree.ptx")
    endif()
    list(APPEND YBI_COMPILE_DEFINITIONS
        YBI_USE_PATH_GUIDING
        YBI_KD_TREE_KERNEL_PATH=\"${YBI_KD_TREE_KERNEL_PATH}\")
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
    list(APPEND YBI_SOURCE
        ${CMAKE_CURRENT_SOURCE_DIR}/third_party/imgui/imgui.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/third_party/imgui/imgui_draw.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/third_party/imgui/imgui_tables.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/third_party/imgui/imgui_widgets.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/third_party/imgui/backends/imgui_impl_glfw.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/third_party/imgui/backends/imgui_impl_vulkan.cpp)
    list(APPEND YBI_THIRD_PARTY_INCLUDE_DIRECTORIES
        "${CMAKE_CURRENT_SOURCE_DIR}/third_party/imgui"
        "${CMAKE_CURRENT_SOURCE_DIR}/third_party/imgui/backends")
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
               ${CMAKE_CURRENT_SOURCE_DIR}/src/material/shader_compilation.cpp
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
target_compile_options(yeoubi_virtual_textures PRIVATE ${YBI_COMPILE_OPTIONS})

# Generate pmj tables exe
add_executable(yeoubi_pmj_tables
                ${CMAKE_CURRENT_SOURCE_DIR}/src/cmd/generate_pmj02_tables.cpp
                ${CMAKE_CURRENT_SOURCE_DIR}/src/render/sampling/pmj02.cpp
)

target_include_directories(yeoubi_pmj_tables PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/src
)

target_compile_definitions(yeoubi_pmj_tables PRIVATE
    ${YBI_COMPILE_DEFINITIONS}
    $<$<CONFIG:Debug,RelWithDebInfo>:YBI_DEBUG>
)
target_compile_options(yeoubi_pmj_tables PRIVATE ${YBI_COMPILE_OPTIONS})

# Generate vMF split table exe
add_executable(generate_vmf_split_table
    ${CMAKE_CURRENT_SOURCE_DIR}/src/cmd/generate_vmf_split_table.cpp
)

target_include_directories(generate_vmf_split_table PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/src
)

target_compile_definitions(generate_vmf_split_table PRIVATE
    ${YBI_COMPILE_DEFINITIONS}
    $<$<CONFIG:Debug,RelWithDebInfo>:YBI_DEBUG>
)
target_compile_options(generate_vmf_split_table PRIVATE ${YBI_COMPILE_OPTIONS})

# Renderer exe
add_executable(yeoubi
               ${CMAKE_CURRENT_SOURCE_DIR}/src/core/allocation.cpp
               ${CMAKE_CURRENT_SOURCE_DIR}/src/io/camera_json.cpp
               ${CMAKE_CURRENT_SOURCE_DIR}/src/io/usd.cpp
               ${CMAKE_CURRENT_SOURCE_DIR}/src/io/usd_material.cpp
               ${CMAKE_CURRENT_SOURCE_DIR}/src/material/shader_compilation.cpp
               ${CMAKE_CURRENT_SOURCE_DIR}/src/material/material_compilation.cpp
               ${CMAKE_CURRENT_SOURCE_DIR}/src/scene/scene.cpp
               ${CMAKE_CURRENT_SOURCE_DIR}/src/render/render_system.cpp
               ${CMAKE_CURRENT_SOURCE_DIR}/src/render/wavefront_integrator.cpp
               ${CMAKE_CURRENT_SOURCE_DIR}/src/render/path_guiding/lut.cpp
               ${CMAKE_CURRENT_SOURCE_DIR}/src/device/device.cpp
               ${CMAKE_CURRENT_SOURCE_DIR}/src/display_device/display_device.cpp
               ${CMAKE_CURRENT_SOURCE_DIR}/src/input/camera_controller.cpp
               ${CMAKE_CURRENT_SOURCE_DIR}/src/render/virtual_texture/manager.cpp
               ${CMAKE_CURRENT_SOURCE_DIR}/src/scene/clusterizer.cpp
               ${CMAKE_CURRENT_SOURCE_DIR}/src/material/lut.cpp
               ${CMAKE_CURRENT_SOURCE_DIR}/src/scene/subdivision.cpp
               ${CMAKE_CURRENT_SOURCE_DIR}/src/device/cpu_device.cpp
               ${CMAKE_CURRENT_SOURCE_DIR}/src/application.cpp
               ${CMAKE_CURRENT_SOURCE_DIR}/src/main.cpp
               ${YBI_SOURCE}
)

if (YBI_USE_PATH_GUIDING)
    target_sources(yeoubi PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/src/kernel/backend/cpu/kd_tree.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/src/render/path_guiding/path_guiding.cpp)
endif()

# Keep the precomputed guiding LUT beside the renderer executable.
set_property(TARGET yeoubi APPEND PROPERTY LINK_DEPENDS
    "${CMAKE_CURRENT_SOURCE_DIR}/assets/lut/guiding/vmf_split_table.bin")
add_custom_command(TARGET yeoubi POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy_if_different
        "${CMAKE_CURRENT_SOURCE_DIR}/assets/lut/guiding/vmf_split_table.bin"
        "$<TARGET_FILE_DIR:yeoubi>/vmf_split_table.bin"
    VERBATIM
)

list(APPEND YBI_THIRD_PARTY_INCLUDE_DIRECTORIES
    "${CMAKE_CURRENT_SOURCE_DIR}/third_party"
)

function (compile_cuda src)
    cmake_parse_arguments(PARSE_ARGV 1 arg
        ""
        "OUTPUT;DEPFILE;ARCHITECTURE"
        "DEFINITIONS;DEPENDS;INCLUDE_DIRS;OPTIONS"
    )

    set(NVCC_OPTIONS ${arg_OPTIONS})
    foreach(definition IN LISTS arg_DEFINITIONS)
        list(APPEND NVCC_OPTIONS "-D${definition}")
    endforeach()
    foreach(include_dir IN LISTS arg_INCLUDE_DIRS)
        list(APPEND NVCC_OPTIONS "-I${include_dir}")
    endforeach()

    add_custom_command(
        OUTPUT "${arg_OUTPUT}"
        COMMAND "${CMAKE_CUDA_COMPILER}"
        --ptx
        "--gpu-architecture=${arg_ARCHITECTURE}"
        "--std=c++17"
        "--compiler-bindir=${YBI_CUDA_HOST_COMPILER_DIR}"
        "--use_fast_math"
        "-I${CMAKE_CURRENT_SOURCE_DIR}/src"
        ${NVCC_OPTIONS}
        -MD
        -MF "${arg_DEPFILE}"
        "-o" "${arg_OUTPUT}"
        "${src}"
        DEPENDS "${src}" ${arg_DEPENDS}
        DEPFILE "${arg_DEPFILE}"
        VERBATIM
    )
endfunction()


if (YBI_HAS_CUDA AND OPTIX_FOUND)
    get_filename_component(YBI_CUDA_HOST_COMPILER_DIR "${CMAKE_CXX_COMPILER}" DIRECTORY)

    set(YBI_CUDA_ARCHITECTURE "compute_75")
    set(YBI_CUDA_STANDARD "c++${CMAKE_CXX_STANDARD}")
    list(GET CUDAToolkit_INCLUDE_DIRS 0 YBI_CUDA_INCLUDE_PATH)
    list(LENGTH CUDAToolkit_INCLUDE_DIRS YBI_CUDA_INCLUDE_PATH_COUNT)
    if (YBI_CUDA_INCLUDE_PATH_COUNT GREATER 1)
        list(GET CUDAToolkit_INCLUDE_DIRS 1 YBI_CUDA_CCCL_INCLUDE_PATH)
    endif()

    set(YBI_OPTIX_OPENPBR_PTX_PATH "${CMAKE_CURRENT_BINARY_DIR}/kernel_openpbr.ptx")
    set(YBI_OPTIX_MATERIALX_PTX_PATH "${CMAKE_CURRENT_BINARY_DIR}/kernel_materialx.ptx")
    set(YBI_OPTIX_INCLUDE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/src")
    set(YBI_OPTIX_SDK_INCLUDE_PATH "${OPTIX_INCLUDE_DIR}")
    foreach(YBI_MATERIAL_BACKEND IN ITEMS OPENPBR MATERIALX)
        string(TOLOWER "${YBI_MATERIAL_BACKEND}" YBI_MATERIAL_BACKEND_LOWER)
        set(YBI_BACKEND_PTX_PATH
            "${CMAKE_CURRENT_BINARY_DIR}/kernel_${YBI_MATERIAL_BACKEND_LOWER}.ptx")
        set(YBI_BACKEND_PTX_DEPFILE "${YBI_BACKEND_PTX_PATH}.d")

        compile_cuda("${CMAKE_CURRENT_SOURCE_DIR}/src/kernel/backend/optix/kernel.cu"
            OUTPUT "${YBI_BACKEND_PTX_PATH}" DEPFILE "${YBI_BACKEND_PTX_DEPFILE}"
            ARCHITECTURE ${YBI_CUDA_ARCHITECTURE}
            DEFINITIONS "OPTIX_VERSION=${OPTIX_VERSION}" "YBI_MATERIAL_BACKEND_${YBI_MATERIAL_BACKEND}=1"
            INCLUDE_DIRS "${YBI_OPTIX_SDK_INCLUDE_PATH}"
            # OPTIONS -lineinfo
            DEPENDS ${YBI_MATERIALX_STUBS}
        )
    endforeach()
    add_custom_target(yeoubi_optix_kernel DEPENDS
        "${YBI_OPTIX_OPENPBR_PTX_PATH}" "${YBI_OPTIX_MATERIALX_PTX_PATH}")
    add_dependencies(yeoubi yeoubi_optix_kernel)

    set(YBI_VIRTUAL_TEXTURE_KERNEL_PATH "${CMAKE_CURRENT_BINARY_DIR}/virtual_textures.ptx")
    set(YBI_VIRTUAL_TEXTURE_KERNEL_DEPFILE
        "${CMAKE_CURRENT_BINARY_DIR}/virtual_textures.ptx.d")
    compile_cuda("${CMAKE_CURRENT_SOURCE_DIR}/src/kernel/backend/cuda/virtual_textures.cu"
        OUTPUT "${YBI_VIRTUAL_TEXTURE_KERNEL_PATH}"
        DEPFILE "${YBI_VIRTUAL_TEXTURE_KERNEL_DEPFILE}"
        ARCHITECTURE sm_75
    )
    add_custom_target(yeoubi_virtual_texture_kernel DEPENDS ${YBI_VIRTUAL_TEXTURE_KERNEL_PATH})
    add_dependencies(yeoubi yeoubi_virtual_texture_kernel)

    # Accumulation
    set(YBI_ACCUMULATE_KERNEL_PATH "${CMAKE_CURRENT_BINARY_DIR}/accumulate.ptx")
    set(YBI_ACCUMULATE_KERNEL_DEPFILE "${CMAKE_CURRENT_BINARY_DIR}/accumulate.ptx.d")

    compile_cuda("${CMAKE_CURRENT_SOURCE_DIR}/src/kernel/backend/cuda/accumulate.cu"
        OUTPUT "${YBI_ACCUMULATE_KERNEL_PATH}"
        DEPFILE "${YBI_ACCUMULATE_KERNEL_DEPFILE}"
        ARCHITECTURE sm_75
    )
    add_custom_target(yeoubi_accumulate_kernel DEPENDS ${YBI_ACCUMULATE_KERNEL_PATH})
    add_dependencies(yeoubi yeoubi_accumulate_kernel)

    # Camera ray generation
    set(YBI_GENERATE_CAMERA_RAY_KERNEL_PATH "${CMAKE_CURRENT_BINARY_DIR}/generate_camera_ray.ptx")
    set(YBI_GENERATE_CAMERA_RAY_KERNEL_DEPFILE "${YBI_GENERATE_CAMERA_RAY_KERNEL_PATH}.d")

    compile_cuda("${CMAKE_CURRENT_SOURCE_DIR}/src/kernel/backend/cuda/generate_camera_ray.cu"
        OUTPUT "${YBI_GENERATE_CAMERA_RAY_KERNEL_PATH}"
        DEPFILE "${YBI_GENERATE_CAMERA_RAY_KERNEL_DEPFILE}"
        ARCHITECTURE sm_75
    )
    add_custom_target(yeoubi_generate_camera_ray_kernel DEPENDS ${YBI_GENERATE_CAMERA_RAY_KERNEL_PATH})
    add_dependencies(yeoubi yeoubi_generate_camera_ray_kernel)

    # Kernel state counts
    set(YBI_COUNT_KERNELS_KERNEL_PATH "${CMAKE_CURRENT_BINARY_DIR}/count_kernels.ptx")
    set(YBI_COUNT_KERNELS_KERNEL_DEPFILE "${YBI_COUNT_KERNELS_KERNEL_PATH}.d")

    compile_cuda("${CMAKE_CURRENT_SOURCE_DIR}/src/kernel/backend/cuda/count_kernels.cu"
        OUTPUT "${YBI_COUNT_KERNELS_KERNEL_PATH}"
        DEPFILE "${YBI_COUNT_KERNELS_KERNEL_DEPFILE}"
        ARCHITECTURE sm_75
    )
    add_custom_target(yeoubi_count_kernels_kernel DEPENDS ${YBI_COUNT_KERNELS_KERNEL_PATH})
    add_dependencies(yeoubi yeoubi_count_kernels_kernel)

    # Sample indices for the selected kernel
    set(YBI_SELECT_SAMPLE_INDICES_KERNEL_PATH "${CMAKE_CURRENT_BINARY_DIR}/select_sample_indices.ptx")
    set(YBI_SELECT_SAMPLE_INDICES_KERNEL_DEPFILE "${YBI_SELECT_SAMPLE_INDICES_KERNEL_PATH}.d")

    compile_cuda("${CMAKE_CURRENT_SOURCE_DIR}/src/kernel/backend/cuda/select_sample_indices.cu"
        OUTPUT "${YBI_SELECT_SAMPLE_INDICES_KERNEL_PATH}"
        DEPFILE "${YBI_SELECT_SAMPLE_INDICES_KERNEL_DEPFILE}"
        ARCHITECTURE sm_75
    )
    add_custom_target(yeoubi_select_sample_indices_kernel DEPENDS ${YBI_SELECT_SAMPLE_INDICES_KERNEL_PATH})
    add_dependencies(yeoubi yeoubi_select_sample_indices_kernel)

    # VMM fitting: one CUDA block per mixture.
    set(YBI_UPDATE_VMM_KERNEL_PATH "${CMAKE_CURRENT_BINARY_DIR}/update_vmm.ptx")
    set(YBI_UPDATE_VMM_KERNEL_DEPFILE "${YBI_UPDATE_VMM_KERNEL_PATH}.d")

    compile_cuda("${CMAKE_CURRENT_SOURCE_DIR}/src/kernel/backend/cuda/update_vmm.cu"
        OUTPUT "${YBI_UPDATE_VMM_KERNEL_PATH}"
        DEPFILE "${YBI_UPDATE_VMM_KERNEL_DEPFILE}"
        ARCHITECTURE sm_75
    )
    add_custom_target(yeoubi_update_vmm_kernel DEPENDS ${YBI_UPDATE_VMM_KERNEL_PATH})
    add_dependencies(yeoubi yeoubi_update_vmm_kernel)

    if (YBI_USE_PATH_GUIDING)
        set(YBI_KD_TREE_KERNEL_DEPFILE "${YBI_KD_TREE_KERNEL_PATH}.d")
        compile_cuda("${CMAKE_CURRENT_SOURCE_DIR}/src/kernel/backend/cuda/kd_tree.cu"
            OUTPUT "${YBI_KD_TREE_KERNEL_PATH}"
            DEPFILE "${YBI_KD_TREE_KERNEL_DEPFILE}"
            ARCHITECTURE sm_75
        )
        add_custom_target(yeoubi_kd_tree_kernel DEPENDS ${YBI_KD_TREE_KERNEL_PATH})
        add_dependencies(yeoubi yeoubi_kd_tree_kernel)
    endif()

    # Stream compaction
    set(YBI_STREAM_COMPACTION_KERNEL_PATH "${CMAKE_CURRENT_BINARY_DIR}/stream_compaction.ptx")
    set(YBI_STREAM_COMPACTION_KERNEL_DEPFILE "${YBI_STREAM_COMPACTION_KERNEL_PATH}.d")

    compile_cuda("${CMAKE_CURRENT_SOURCE_DIR}/src/kernel/backend/cuda/stream_compaction.cu"
        OUTPUT "${YBI_STREAM_COMPACTION_KERNEL_PATH}"
        DEPFILE "${YBI_STREAM_COMPACTION_KERNEL_DEPFILE}"
        ARCHITECTURE sm_75
    )
    add_custom_target(yeoubi_stream_compaction_kernel DEPENDS ${YBI_STREAM_COMPACTION_KERNEL_PATH})
    add_dependencies(yeoubi yeoubi_stream_compaction_kernel)

    # Wavefront OptiX kernels
    set(YBI_WAVEFRONT_KERNEL_PATH "${CMAKE_CURRENT_BINARY_DIR}/kernels.ptx")
    set(YBI_WAVEFRONT_KERNEL_DEPFILE "${YBI_WAVEFRONT_KERNEL_PATH}.d")

    compile_cuda("${CMAKE_CURRENT_SOURCE_DIR}/src/kernel/backend/optix/kernels.cu"
        OUTPUT "${YBI_WAVEFRONT_KERNEL_PATH}"
        DEPFILE "${YBI_WAVEFRONT_KERNEL_DEPFILE}"
        ARCHITECTURE ${YBI_CUDA_ARCHITECTURE}
        DEFINITIONS "OPTIX_VERSION=${OPTIX_VERSION}"
        INCLUDE_DIRS "${YBI_OPTIX_SDK_INCLUDE_PATH}"
    )
    add_custom_target(yeoubi_wavefront_kernel DEPENDS ${YBI_WAVEFRONT_KERNEL_PATH})
    add_dependencies(yeoubi yeoubi_wavefront_kernel)

    list (APPEND YBI_COMPILE_DEFINITIONS
        YBI_WAVEFRONT_KERNEL_PATH=\"${YBI_WAVEFRONT_KERNEL_PATH}\"
        YBI_OPTIX_OPENPBR_PTX_PATH=\"${YBI_OPTIX_OPENPBR_PTX_PATH}\"
        YBI_OPTIX_MATERIALX_PTX_PATH=\"${YBI_OPTIX_MATERIALX_PTX_PATH}\")
    list (APPEND YBI_COMPILE_DEFINITIONS YBI_VIRTUAL_TEXTURE_KERNEL_PATH=\"${YBI_VIRTUAL_TEXTURE_KERNEL_PATH}\")
    list (APPEND YBI_COMPILE_DEFINITIONS YBI_ACCUMULATE_KERNEL_PATH=\"${YBI_ACCUMULATE_KERNEL_PATH}\")
    list (APPEND YBI_COMPILE_DEFINITIONS YBI_GENERATE_CAMERA_RAY_KERNEL_PATH=\"${YBI_GENERATE_CAMERA_RAY_KERNEL_PATH}\")
    list (APPEND YBI_COMPILE_DEFINITIONS YBI_COUNT_KERNELS_KERNEL_PATH=\"${YBI_COUNT_KERNELS_KERNEL_PATH}\")
    list (APPEND YBI_COMPILE_DEFINITIONS YBI_SELECT_SAMPLE_INDICES_KERNEL_PATH=\"${YBI_SELECT_SAMPLE_INDICES_KERNEL_PATH}\")
    list (APPEND YBI_COMPILE_DEFINITIONS YBI_UPDATE_VMM_KERNEL_PATH=\"${YBI_UPDATE_VMM_KERNEL_PATH}\")
    list (APPEND YBI_COMPILE_DEFINITIONS YBI_STREAM_COMPACTION_KERNEL_PATH=\"${YBI_STREAM_COMPACTION_KERNEL_PATH}\")
    list (APPEND YBI_COMPILE_DEFINITIONS YBI_OPTIX_INCLUDE_PATH=\"${YBI_OPTIX_INCLUDE_PATH}\")
    list (APPEND YBI_COMPILE_DEFINITIONS YBI_OPTIX_SDK_INCLUDE_PATH=\"${YBI_OPTIX_SDK_INCLUDE_PATH}\")
    list (APPEND YBI_COMPILE_DEFINITIONS YBI_CUDA_INCLUDE_PATH=\"${YBI_CUDA_INCLUDE_PATH}\")
    if (YBI_CUDA_CCCL_INCLUDE_PATH)
        list (APPEND YBI_COMPILE_DEFINITIONS YBI_CUDA_CCCL_INCLUDE_PATH=\"${YBI_CUDA_CCCL_INCLUDE_PATH}\")
    endif()
    list (APPEND YBI_COMPILE_DEFINITIONS YBI_CUDA_ARCHITECTURE=\"${YBI_CUDA_ARCHITECTURE}\")
    list (APPEND YBI_COMPILE_DEFINITIONS YBI_CUDA_STANDARD=\"${YBI_CUDA_STANDARD}\")
    list (APPEND YBI_COMPILE_DEFINITIONS YBI_OPTIX_VERSION_STRING=\"${OPTIX_VERSION}\")
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
target_compile_options(yeoubi PRIVATE ${YBI_COMPILE_OPTIONS})

# add_executable(material_compilation_dielectric_test
#     ${CMAKE_CURRENT_SOURCE_DIR}/tests/material_compilation_dielectric_test.cpp
#     ${CMAKE_CURRENT_SOURCE_DIR}/src/core/allocation.cpp
#     ${CMAKE_CURRENT_SOURCE_DIR}/src/device/device.cpp
#     ${CMAKE_CURRENT_SOURCE_DIR}/src/display_device/display_device.cpp
#     ${CMAKE_CURRENT_SOURCE_DIR}/src/device/cuda_device.cpp
#     ${CMAKE_CURRENT_SOURCE_DIR}/src/kernel/backend/cuda/clusters.cu
#     ${CMAKE_CURRENT_SOURCE_DIR}/src/material/material_compilation.cpp
#     ${CMAKE_CURRENT_SOURCE_DIR}/src/scene/clusterizer.cpp
#     ${CMAKE_CURRENT_SOURCE_DIR}/src/scene/scene.cpp
# )
#
# target_include_directories(material_compilation_dielectric_test PRIVATE
#     ${CMAKE_CURRENT_SOURCE_DIR}/src
#     ${OPENUSD_INCLUDE_DIRS}
#     ${YBI_INCLUDE_DIRS}
# )
#
# target_include_directories(material_compilation_dielectric_test SYSTEM PRIVATE
#     ${YBI_THIRD_PARTY_INCLUDE_DIRECTORIES}
# )
#
# target_link_libraries(material_compilation_dielectric_test PRIVATE
#     ${YBI_DEPENDENCIES}
# )
#
# target_compile_definitions(material_compilation_dielectric_test PRIVATE
#     ${YBI_COMPILE_DEFINITIONS}
#     $<$<CONFIG:Debug,RelWithDebInfo>:YBI_DEBUG>
# )
