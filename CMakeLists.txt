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
    endif()
endif()
# Hardware Specific End

add_executable(yeoubi ${YBI_SOURCE})

target_compile_definitions(yeoubi PRIVATE ${YBI_COMPILE_DEFINITIONS})
target_include_directories(yeoubi PRIVATE ${YBI_INCLUDE_DIRECTORIES})
target_include_directories(yeoubi SYSTEM PRIVATE ${YBI_THIRD_PARTY_INCLUDE_DIRECTORIES})
target_link_libraries(yeoubi PRIVATE ${YBI_DEPENDENCIES})
