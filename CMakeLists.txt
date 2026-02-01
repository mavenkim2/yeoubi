cmake_minimum_required(VERSION 3.10)
project (yeoubi)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/src/cmake")

set (YBI_DEPENDENCIES "")
set (YBI_INCLUDE_DIRECTORIES "")
set (YBI_COMPILE_DEFINITIONS 
    "YBI_NAMESPACE_BEGIN=namespace ybi {"
    "YBI_NAMESPACE_END=}"
)

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
        list (APPEND YBI_INCLUDE_DIRECTORIES "${OPTIX_INCLUDE_DIR}")
    endif()
endif()

add_executable(yeoubi ${CMAKE_CURRENT_SOURCE_DIR}/src/test.cpp)

target_link_libraries(yeoubi PRIVATE ${YBI_DEPENDENCIES})
target_compile_definitions(yeoubi PRIVATE ${YBI_COMPILE_DEFINITIONS})
target_include_directories(yeoubi PRIVATE ${YBI_INCLUDE_DIRECTORIES})
