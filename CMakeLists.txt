cmake_minimum_required(VERSION 3.10)
project (yeoubi)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
option(WITH_EMBREE "Enable Embree CPU backend" ON)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/src/cmake")
include(CheckLanguage)
check_language(CUDA)
set(YBI_HAS_CUDA OFF)
if (CMAKE_CUDA_COMPILER)
    set(YBI_HAS_CUDA ON)
endif()

set (YBI_DEPENDENCIES "")
set (YBI_INCLUDE_DIRECTORIES "")
set (YBI_THIRD_PARTY_INCLUDE_DIRECTORIES "")
set (YBI_COMPILE_DEFINITIONS 
    "YBI_NAMESPACE_BEGIN=namespace ybi {"
    "YBI_NAMESPACE_END=}"
)
set (YBI_SOURCE "")

list (APPEND YBI_INCLUDE_DIRECTORIES ${CMAKE_CURRENT_SOURCE_DIR}/src)
list (APPEND YBI_SOURCE 
    ${CMAKE_CURRENT_SOURCE_DIR}/src/test.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/device/device.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/scene/scene.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/io/usd/load.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/io/usd/instance_dag_build.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/scene/attributes.cpp
)
list (APPEND YBI_COMPILE_DEFINITIONS 
    $<$<CONFIG:Debug,RelWithDebInfo>:YBI_DEBUG>)


# Compiler specific
if (MSVC)
    add_compile_definitions(_USE_MATH_DEFINES)
endif()

# Add third party libraries
add_subdirectory(src/third_party)

if (WITH_USD)
    list (APPEND YBI_DEPENDENCIES ${OPENUSD_LIBRARIES})
    list (APPEND YBI_THIRD_PARTY_INCLUDE_DIRECTORIES ${OPENUSD_INCLUDE_DIRS})
    # TODO IMPORTANT: right now on my local machine I'm just manually setting the 
    # PATH for the openusd download. handle this for other machines
    # install(FILES ${OPENUSD_LIBRARIES} DESTINATION ${CMAKE_INSTALL_PREFIX})
endif()

if (WITH_USD)
    add_executable(yeoubi_usd_add_look_camera
        ${CMAKE_CURRENT_SOURCE_DIR}/cmd/usd_add_look_camera.cpp)
    target_compile_definitions(yeoubi_usd_add_look_camera PRIVATE ${YBI_COMPILE_DEFINITIONS})
    if (WIN32)
        target_compile_definitions(yeoubi_usd_add_look_camera PRIVATE NOMINMAX WIN32_LEAN_AND_MEAN)
    endif()
    target_include_directories(yeoubi_usd_add_look_camera PRIVATE ${YBI_INCLUDE_DIRECTORIES})
    target_include_directories(yeoubi_usd_add_look_camera SYSTEM PRIVATE ${YBI_THIRD_PARTY_INCLUDE_DIRECTORIES})
    target_link_libraries(yeoubi_usd_add_look_camera PRIVATE ${OPENUSD_LIBRARIES})

    add_executable(yeoubi_usd_extract
        ${CMAKE_CURRENT_SOURCE_DIR}/tests/bvh/usd_extract.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/tests/bvh/usd_subdiv_select.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/tests/io/usd_mesh_io.cpp)
    target_compile_definitions(yeoubi_usd_extract PRIVATE ${YBI_COMPILE_DEFINITIONS})
    if (WIN32)
        target_compile_definitions(yeoubi_usd_extract PRIVATE NOMINMAX WIN32_LEAN_AND_MEAN)
    endif()
    target_include_directories(yeoubi_usd_extract PRIVATE
        ${YBI_INCLUDE_DIRECTORIES}
        ${CMAKE_CURRENT_SOURCE_DIR}/tests)
    target_include_directories(yeoubi_usd_extract SYSTEM PRIVATE ${YBI_THIRD_PARTY_INCLUDE_DIRECTORIES})
    target_link_libraries(yeoubi_usd_extract PRIVATE ${OPENUSD_LIBRARIES})

    add_executable(yeoubi_tessellation_adaptive_impl
        ${CMAKE_CURRENT_SOURCE_DIR}/tests/tesselation/tessellation_adaptive_impl.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/src/tessellation/subdivision.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/src/tessellation/edge_map_validation.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/src/io/usd/load.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/src/io/usd/instance_dag_build.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/src/scene/scene.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/src/scene/attributes.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/tests/bvh/usd_subdiv_select.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/tests/io/usd_subdiv_json_io.cpp)
    target_compile_definitions(yeoubi_tessellation_adaptive_impl PRIVATE ${YBI_COMPILE_DEFINITIONS})
    if (WIN32)
        target_compile_definitions(yeoubi_tessellation_adaptive_impl PRIVATE NOMINMAX WIN32_LEAN_AND_MEAN)
    endif()
    target_include_directories(yeoubi_tessellation_adaptive_impl PRIVATE
        ${YBI_INCLUDE_DIRECTORIES}
        ${CMAKE_CURRENT_SOURCE_DIR}/tests)
    target_include_directories(yeoubi_tessellation_adaptive_impl SYSTEM PRIVATE ${YBI_THIRD_PARTY_INCLUDE_DIRECTORIES})
    target_link_libraries(yeoubi_tessellation_adaptive_impl PRIVATE ${OPENUSD_LIBRARIES})

    add_executable(yeoubi_usd_instance_dag_test
        ${CMAKE_CURRENT_SOURCE_DIR}/tests/io/usd_instance_dag_test.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/src/io/usd/instance_dag_build.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/src/io/usd/load.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/src/scene/scene.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/src/scene/attributes.cpp)
    target_compile_definitions(yeoubi_usd_instance_dag_test PRIVATE ${YBI_COMPILE_DEFINITIONS})
    if (WIN32)
        target_compile_definitions(yeoubi_usd_instance_dag_test PRIVATE NOMINMAX WIN32_LEAN_AND_MEAN)
    endif()
    target_include_directories(yeoubi_usd_instance_dag_test PRIVATE
        ${YBI_INCLUDE_DIRECTORIES}
        ${CMAKE_CURRENT_SOURCE_DIR}/tests)
    target_include_directories(yeoubi_usd_instance_dag_test SYSTEM PRIVATE ${YBI_THIRD_PARTY_INCLUDE_DIRECTORIES})
    target_link_libraries(yeoubi_usd_instance_dag_test PRIVATE ${OPENUSD_LIBRARIES})

    add_executable(yeoubi_usd_tile_prepare
        ${CMAKE_CURRENT_SOURCE_DIR}/cmd/virtual_texture/main.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/cmd/usd_ntc_encode/cli.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/cmd/usd_ntc_encode/usd.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/cmd/virtual_texture/tiles.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/cmd/virtual_texture/tile_binary.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/src/texture/exr_io.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/src/texture/path_utils.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/src/texture/udim_utils.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/src/third_party/tinyexr/deps/miniz/miniz.c)
    target_compile_definitions(yeoubi_usd_tile_prepare PRIVATE ${YBI_COMPILE_DEFINITIONS})
    if (WIN32)
        target_compile_definitions(yeoubi_usd_tile_prepare PRIVATE NOMINMAX WIN32_LEAN_AND_MEAN)
    endif()
    target_include_directories(yeoubi_usd_tile_prepare PRIVATE
        ${YBI_INCLUDE_DIRECTORIES}
        ${CMAKE_CURRENT_SOURCE_DIR}/src/third_party/tinyexr
        ${CMAKE_CURRENT_SOURCE_DIR}/src/third_party/tinyexr/deps/miniz
        ${CMAKE_CURRENT_SOURCE_DIR}/tests)
    target_include_directories(yeoubi_usd_tile_prepare SYSTEM PRIVATE ${YBI_THIRD_PARTY_INCLUDE_DIRECTORIES})
    target_link_libraries(yeoubi_usd_tile_prepare PRIVATE ${OPENUSD_LIBRARIES})

    if (YBI_HAS_CUDA AND WITH_NTC)
        add_executable(yeoubi_usd_ntc_encode
            ${CMAKE_CURRENT_SOURCE_DIR}/cmd/usd_ntc_encode/main.cpp
            ${CMAKE_CURRENT_SOURCE_DIR}/cmd/usd_ntc_encode/cli.cpp
            ${CMAKE_CURRENT_SOURCE_DIR}/cmd/usd_ntc_encode/usd.cpp
            ${CMAKE_CURRENT_SOURCE_DIR}/cmd/usd_ntc_encode/ntc.cpp
            ${CMAKE_CURRENT_SOURCE_DIR}/src/texture/exr_io.cpp
            ${CMAKE_CURRENT_SOURCE_DIR}/src/texture/path_utils.cpp
            ${CMAKE_CURRENT_SOURCE_DIR}/src/texture/udim_utils.cpp
            ${CMAKE_CURRENT_SOURCE_DIR}/cmd/virtual_texture/tiles.cpp
            ${CMAKE_CURRENT_SOURCE_DIR}/cmd/virtual_texture/tile_binary.cpp
            ${CMAKE_CURRENT_SOURCE_DIR}/src/third_party/tinyexr/deps/miniz/miniz.c)
        target_compile_definitions(yeoubi_usd_ntc_encode PRIVATE ${YBI_COMPILE_DEFINITIONS})
        if (WIN32)
            target_compile_definitions(yeoubi_usd_ntc_encode PRIVATE NOMINMAX WIN32_LEAN_AND_MEAN)
        endif()
        target_include_directories(yeoubi_usd_ntc_encode PRIVATE
            ${YBI_INCLUDE_DIRECTORIES}
            ${NTC_INCLUDE_DIRS}
            ${CMAKE_CURRENT_SOURCE_DIR}/tests)
        target_include_directories(yeoubi_usd_ntc_encode SYSTEM PRIVATE ${YBI_THIRD_PARTY_INCLUDE_DIRECTORIES})
        target_link_libraries(yeoubi_usd_ntc_encode PRIVATE ${OPENUSD_LIBRARIES} ${NTC_LIBRARIES})

        add_executable(yeoubi_tile_verify
            ${CMAKE_CURRENT_SOURCE_DIR}/cmd/virtual_texture/tile_verify.cpp
            ${CMAKE_CURRENT_SOURCE_DIR}/cmd/virtual_texture/tile_binary.cpp
            ${CMAKE_CURRENT_SOURCE_DIR}/src/texture/exr_io.cpp
            ${CMAKE_CURRENT_SOURCE_DIR}/src/texture/udim_utils.cpp
            ${CMAKE_CURRENT_SOURCE_DIR}/src/third_party/tinyexr/deps/miniz/miniz.c)
        target_compile_definitions(yeoubi_tile_verify PRIVATE ${YBI_COMPILE_DEFINITIONS})
        if (WIN32)
            target_compile_definitions(yeoubi_tile_verify PRIVATE NOMINMAX WIN32_LEAN_AND_MEAN)
        endif()
        target_include_directories(yeoubi_tile_verify PRIVATE
            ${YBI_INCLUDE_DIRECTORIES}
            ${CMAKE_CURRENT_SOURCE_DIR}/src/third_party/tinyexr
            ${CMAKE_CURRENT_SOURCE_DIR}/src/third_party/tinyexr/deps/miniz)
    else()
        message(STATUS "CUDA/NTC not found; yeoubi_usd_ntc_encode disabled.")
    endif()
endif()

if (MESHOPTIMIZER_FOUND)
    list (APPEND YBI_DEPENDENCIES ${MESHOPTIMIZER_LIBRARY})
    list (APPEND YBI_THIRD_PARTY_INCLUDE_DIRECTORIES ${MESHOPTIMIZER_INCLUDE_DIR})
endif()

add_executable(yeoubi_render
    ${CMAKE_CURRENT_SOURCE_DIR}/cmd/render/main.cpp)

if (WITH_EMBREE)
    find_package(Embree MODULE)
    if (EMBREE_FOUND)
        list (APPEND YBI_DEPENDENCIES ${EMBREE_LIBRARIES})
        list (APPEND YBI_THIRD_PARTY_INCLUDE_DIRECTORIES ${EMBREE_INCLUDE_DIRS})
        list (APPEND YBI_COMPILE_DEFINITIONS WITH_EMBREE YBI_EMBREE_VERSION_MAJOR=${EMBREE_VERSION_MAJOR})
        list (APPEND YBI_SOURCE ${CMAKE_CURRENT_SOURCE_DIR}/src/device/cpu_device.cpp)
        list (APPEND YBI_SOURCE ${CMAKE_CURRENT_SOURCE_DIR}/src/device/bvh/embree/bvh.cpp)
        list (APPEND YBI_SOURCE ${CMAKE_CURRENT_SOURCE_DIR}/src/render/cpu_kernels.cpp)
        message (STATUS "Embree Found: ${EMBREE_INCLUDE_DIRS} ${EMBREE_VERSION}")
    else()
        message(STATUS "Embree not found; CPU Embree backend disabled.")
    endif()
endif()

# OS Specific
if (WIN32)
    list (APPEND YBI_COMPILE_DEFINITIONS "NOMINMAX" "WIN32_LEAN_AND_MEAN")
endif()
# OS Specific End

# Hardware Specific
# TODO: check_language for CUDA only works on MSVC?
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
        list (APPEND YBI_SOURCE ${CMAKE_CURRENT_SOURCE_DIR}/src/device/bvh/optix/bvh.cpp)
        message (STATUS "OptiX Found: ${OPTIX_INCLUDE_DIR} ${OPTIX_VERSION}")
        if (OPTIX_VERSION VERSION_GREATER_EQUAL 90000)
            list (APPEND YBI_SOURCE ${CMAKE_CURRENT_SOURCE_DIR}/src/device/bvh/optix/clusters.cu)
            message (STATUS "OptiX CLAS supported.")
        endif()
    endif()
endif()
# Hardware Specific End

add_executable(yeoubi ${YBI_SOURCE})

target_compile_definitions(yeoubi PRIVATE ${YBI_COMPILE_DEFINITIONS})
target_include_directories(yeoubi PRIVATE ${YBI_INCLUDE_DIRECTORIES})
target_include_directories(yeoubi SYSTEM PRIVATE ${YBI_THIRD_PARTY_INCLUDE_DIRECTORIES})
target_link_libraries(yeoubi PRIVATE ${YBI_DEPENDENCIES})

if (CMAKE_CUDA_COMPILER AND OPTIX_FOUND)
    list(GET CUDAToolkit_INCLUDE_DIRS 0 YBI_CUDA_INCLUDE_DIR)
    set(YBI_OPTIX_PTX_GPU_ARCH "compute_89")
    message(STATUS "OptiX PTX compile arch: ${YBI_OPTIX_PTX_GPU_ARCH}")
    set(YBI_OPTIX_PRIMARY_PTX_PATH "${CMAKE_CURRENT_BINARY_DIR}/optix_primary_programs.ptx")
    set(YBI_OPTIX_PRIMARY_PTX_FROZEN_PATH
        "${CMAKE_CURRENT_SOURCE_DIR}/tests/bvh/optix/frozen/optix_primary_programs_2026_03_03.ptx")
    add_custom_command(
        OUTPUT "${YBI_OPTIX_PRIMARY_PTX_PATH}"
        COMMAND "${CMAKE_CUDA_COMPILER}"
            --ptx
            "--gpu-architecture=${YBI_OPTIX_PTX_GPU_ARCH}"
            "--std=c++14"
            "--use_fast_math"
            "-I${CMAKE_CURRENT_SOURCE_DIR}/src"
            "-I${OPTIX_INCLUDE_DIR}"
            "-I${YBI_CUDA_INCLUDE_DIR}"
            "-I${CMAKE_CURRENT_SOURCE_DIR}/src/third_party/RTXNTC/libraries/RTXNTC-Library/include"
            "-o" "${YBI_OPTIX_PRIMARY_PTX_PATH}"
            "${CMAKE_CURRENT_SOURCE_DIR}/tests/bvh/optix/optix_primary_programs.cu"
        DEPENDS
            "${CMAKE_CURRENT_SOURCE_DIR}/tests/bvh/optix/optix_primary_programs.cu"
        VERBATIM)
    add_custom_target(yeoubi_optix_primary_programs_ptx DEPENDS "${YBI_OPTIX_PRIMARY_PTX_PATH}")

endif()

if (WITH_USD)
    add_executable(yeoubi_template_harness
        ${CMAKE_CURRENT_SOURCE_DIR}/src/render/template_harness.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/src/device/device.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/src/scene/clusterizer.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/src/scene/scene.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/src/texture/path_utils.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/src/texture/exr_io.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/src/texture/udim_utils.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/src/texture/virtual_texture/helpers.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/src/texture/virtual_texture/feedback.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/src/texture/virtual_texture/tile_file.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/src/io/usd/load.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/src/io/usd/instance_dag_build.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/src/scene/attributes.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/src/tessellation/subdivision.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/src/tessellation/edge_map_validation.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/src/third_party/tinyexr/deps/miniz/miniz.c)

    if (WITH_EMBREE AND EMBREE_FOUND)
        target_sources(yeoubi_template_harness PRIVATE
            ${CMAKE_CURRENT_SOURCE_DIR}/src/device/cpu_device.cpp
            ${CMAKE_CURRENT_SOURCE_DIR}/src/device/bvh/embree/bvh.cpp
            ${CMAKE_CURRENT_SOURCE_DIR}/src/render/cpu_kernels.cpp)
    endif()

    if (CMAKE_CUDA_COMPILER AND OPTIX_FOUND)
        target_sources(yeoubi_template_harness PRIVATE
            ${CMAKE_CURRENT_SOURCE_DIR}/src/device/cuda_device.cpp
            ${CMAKE_CURRENT_SOURCE_DIR}/src/device/bvh/optix/bvh.cpp)
        if (OPTIX_VERSION VERSION_GREATER_EQUAL 90000)
            target_sources(yeoubi_template_harness PRIVATE
                ${CMAKE_CURRENT_SOURCE_DIR}/src/device/bvh/optix/clusters.cu)
        endif()

        target_compile_definitions(yeoubi_template_harness PRIVATE
            YBI_OPTIX_PRIMARY_PTX_PATH="${YBI_OPTIX_PRIMARY_PTX_PATH}")
        add_dependencies(yeoubi_template_harness yeoubi_optix_primary_programs_ptx)
    endif()

    target_compile_definitions(yeoubi_template_harness PRIVATE ${YBI_COMPILE_DEFINITIONS})
    target_include_directories(yeoubi_template_harness PRIVATE ${YBI_INCLUDE_DIRECTORIES})
    target_include_directories(yeoubi_template_harness PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/src/third_party/tinyexr
        ${CMAKE_CURRENT_SOURCE_DIR}/src/third_party/tinyexr/deps/miniz)
    target_include_directories(yeoubi_template_harness SYSTEM PRIVATE ${YBI_THIRD_PARTY_INCLUDE_DIRECTORIES})
    target_link_libraries(yeoubi_template_harness PRIVATE ${YBI_DEPENDENCIES})

    if (WITH_NTC AND CMAKE_CUDA_COMPILER AND OPTIX_FOUND)
        target_sources(yeoubi_template_harness PRIVATE
            ${CMAKE_CURRENT_SOURCE_DIR}/tests/bvh/optix/optix_ntc_runtime.cpp)
        target_compile_definitions(yeoubi_template_harness PRIVATE YBI_OPTIX_HARNESS_WITH_NTC=1)
        target_include_directories(yeoubi_template_harness PRIVATE ${NTC_INCLUDE_DIRS})
        target_link_libraries(yeoubi_template_harness PRIVATE ${NTC_LIBRARIES})
    endif()
endif()
