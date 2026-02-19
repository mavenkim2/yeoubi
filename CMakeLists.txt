cmake_minimum_required(VERSION 3.10)
project (yeoubi)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/src/cmake")

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
    ${CMAKE_CURRENT_SOURCE_DIR}/src/scene/scene.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/scene/subdivision.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/io/usd/load.cpp
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

    add_executable(yeoubi_tessellation
        ${CMAKE_CURRENT_SOURCE_DIR}/tests/tesselation/tessellation_adaptive_impl.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/tests/tesselation/edge_rate_obj_io.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/tests/bvh/usd_subdiv_select.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/tests/io/usd_mesh_io.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/tests/io/usd_subdiv_json_io.cpp)
    target_compile_definitions(yeoubi_tessellation PRIVATE ${YBI_COMPILE_DEFINITIONS})
    if (WIN32)
        target_compile_definitions(yeoubi_tessellation PRIVATE NOMINMAX WIN32_LEAN_AND_MEAN)
    endif()
    target_include_directories(yeoubi_tessellation PRIVATE
        ${YBI_INCLUDE_DIRECTORIES}
        ${CMAKE_CURRENT_SOURCE_DIR}/tests)
    target_include_directories(yeoubi_tessellation SYSTEM PRIVATE ${YBI_THIRD_PARTY_INCLUDE_DIRECTORIES})
    target_link_libraries(yeoubi_tessellation PRIVATE ${OPENUSD_LIBRARIES})
endif()

if (MESHOPTIMIZER_FOUND)
    list (APPEND YBI_DEPENDENCIES ${MESHOPTIMIZER_LIBRARY})
    list (APPEND YBI_THIRD_PARTY_INCLUDE_DIRECTORIES ${MESHOPTIMIZER_INCLUDE_DIR})
endif()

# OS Specific
if (WIN32)
    list (APPEND YBI_COMPILE_DEFINITIONS "NOMINMAX" "WIN32_LEAN_AND_MEAN")
endif()
# OS Specific End

# Hardware Specific
# TODO: check_language for CUDA only works on MSVC?
include(CheckLanguage)
check_language(CUDA)
if (CMAKE_CUDA_COMPILER)
    find_package(CUDAToolkit REQUIRED)
    enable_language(CUDA)
    message(STATUS "CUDA Found: ${CMAKE_CUDA_COMPILER} (${CUDA_VERSION})")
    list (APPEND YBI_DEPENDENCIES CUDA::cuda_driver CUDA::cudart)
    list (APPEND YBI_COMPILE_DEFINITIONS WITH_CUDA)

    find_package(OptiX)
    if (OPTIX_FOUND)
        list (APPEND YBI_COMPILE_DEFINITIONS WITH_OPTIX OPTIX_VERSION=${OPTIX_VERSION})
        list (APPEND YBI_THIRD_PARTY_INCLUDE_DIRECTORIES "${OPTIX_INCLUDE_DIR}")
        list (APPEND YBI_SOURCE ${CMAKE_CURRENT_SOURCE_DIR}/src/device/cuda_device.cpp)
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
    set(YBI_OPTIX_PRIMARY_PTX_PATH "${CMAKE_CURRENT_BINARY_DIR}/optix_primary_programs.ptx")
    add_custom_command(
        OUTPUT "${YBI_OPTIX_PRIMARY_PTX_PATH}"
        COMMAND "${CMAKE_CUDA_COMPILER}"
            --ptx
            "--gpu-architecture=compute_75"
            "--std=c++14"
            "--use_fast_math"
            "-I${CMAKE_CURRENT_SOURCE_DIR}/src"
            "-I${OPTIX_INCLUDE_DIR}"
            "-I${YBI_CUDA_INCLUDE_DIR}"
            "-o" "${YBI_OPTIX_PRIMARY_PTX_PATH}"
            "${CMAKE_CURRENT_SOURCE_DIR}/tests/bvh/optix_primary_programs.cu"
        DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/tests/bvh/optix_primary_programs.cu"
        VERBATIM)
    add_custom_target(yeoubi_optix_primary_programs_ptx DEPENDS "${YBI_OPTIX_PRIMARY_PTX_PATH}")

    add_executable(yeoubi_optix_harness
        ${CMAKE_CURRENT_SOURCE_DIR}/tests/bvh/optix_primary_harness.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/src/device/cuda_device.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/src/scene/clusterizer.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/src/scene/scene.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/src/io/usd/load.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/src/scene/attributes.cpp)
    if (OPTIX_VERSION VERSION_GREATER_EQUAL 90000)
        target_sources(yeoubi_optix_harness PRIVATE
            ${CMAKE_CURRENT_SOURCE_DIR}/src/device/bvh/optix/clusters.cu)
    endif()
    add_dependencies(yeoubi_optix_harness yeoubi_optix_primary_programs_ptx)

    target_compile_definitions(yeoubi_optix_harness PRIVATE
        ${YBI_COMPILE_DEFINITIONS}
        YBI_OPTIX_PRIMARY_PTX_PATH="${YBI_OPTIX_PRIMARY_PTX_PATH}")
    target_include_directories(yeoubi_optix_harness PRIVATE ${YBI_INCLUDE_DIRECTORIES})
    target_include_directories(yeoubi_optix_harness SYSTEM PRIVATE ${YBI_THIRD_PARTY_INCLUDE_DIRECTORIES})
    target_link_libraries(yeoubi_optix_harness PRIVATE ${YBI_DEPENDENCIES})
endif()
