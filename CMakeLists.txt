cmake_minimum_required(VERSION 3.16)

project(Yeoubi LANGUAGES CXX) 

set(CMAKE_CXX_STANDARD 17 CACHE STRING "C++ standard for yeoubi and subprojects" FORCE)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/src/cmake")

# Compile Options
option(WITH_CUDA "Enable CUDA backend" OFF)

set(YBI_COMPILE_DEFINITIONS)
set(YBI_DEPENDENCIES)
set(YBI_SOURCE)
if (WIN32)
    list(APPEND YBI_COMPILE_DEFINITIONS NOMINMAX)
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

# Create exe
add_executable(yeoubi
               ${CMAKE_CURRENT_SOURCE_DIR}/src/core/allocation.cpp
               ${CMAKE_CURRENT_SOURCE_DIR}/src/io/usd.cpp
               ${CMAKE_CURRENT_SOURCE_DIR}/src/scene/scene.cpp
               ${CMAKE_CURRENT_SOURCE_DIR}/src/device/device.cpp
               ${CMAKE_CURRENT_SOURCE_DIR}/src/main.cpp
               ${YBI_SOURCE}
)

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
