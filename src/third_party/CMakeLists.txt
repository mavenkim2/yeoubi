# NOTE: this was build by running the provided build script in OpenUSD with --no-python,
# then setting USD_ROOT to the install location.
set(USD_ROOT "" CACHE PATH "Path to the OpenUSD installation directory")
option(WITH_NTC "Enable RTX NTC library" ON)
set(OPENPBR_BSDF_INCLUDE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/openpbr-bsdf")

if (USD_ROOT)
    find_package(OpenUSD)
    if (OPENUSD_FOUND)
        set (WITH_USD ON PARENT_SCOPE) 
        set (OPENUSD_INCLUDE_DIRS ${OPENUSD_INCLUDE_DIRS} PARENT_SCOPE)
        set (OPENUSD_LIBRARIES ${OPENUSD_LIBRARIES} PARENT_SCOPE)
        set (OPENUSD_LIBRARY_DIR ${OPENUSD_LIBRARY_DIR} PARENT_SCOPE)
    endif()
endif()

if (EXISTS "${OPENPBR_BSDF_INCLUDE_DIR}/openpbr.h")
    set(OPENPBR_BSDF_FOUND ON PARENT_SCOPE)
    set(OPENPBR_BSDF_INCLUDE_DIR "${OPENPBR_BSDF_INCLUDE_DIR}" PARENT_SCOPE)
endif()

add_subdirectory(meshoptimizer)

if (TARGET meshoptimizer)
    set (MESHOPTIMIZER_FOUND ON PARENT_SCOPE)
    set(MESHOPTIMIZER_INCLUDE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/meshoptimizer/src" PARENT_SCOPE)
    set(MESHOPTIMIZER_LIBRARY meshoptimizer PARENT_SCOPE)
    message(STATUS "Found meshoptimizer target!")
endif()

if (WITH_NTC AND CMAKE_CUDA_COMPILER
        AND EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/RTXNTC/libraries/RTXNTC-Library/CMakeLists.txt")
    set(NTC_BUILD_SHARED OFF CACHE BOOL "" FORCE)
    set(NTC_WITH_CUDA ON CACHE BOOL "" FORCE)
    set(NTC_WITH_DX12 OFF CACHE BOOL "" FORCE)
    set(NTC_WITH_VULKAN OFF CACHE BOOL "" FORCE)
    set(NTC_WITH_PREBUILT_SHADERS OFF CACHE BOOL "" FORCE)
    add_subdirectory(RTXNTC/libraries/RTXNTC-Library EXCLUDE_FROM_ALL)

    if (TARGET libntc)
        if (CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
            target_compile_options(libntc PRIVATE $<$<COMPILE_LANGUAGE:CXX>:-include cstring>)
        endif()
        set(WITH_NTC ON PARENT_SCOPE)
        set(NTC_LIBRARIES libntc PARENT_SCOPE)
        set(NTC_INCLUDE_DIRS
            "${CMAKE_CURRENT_SOURCE_DIR}/RTXNTC/libraries/RTXNTC-Library/include"
            "${CMAKE_CURRENT_SOURCE_DIR}/tinyexr"
            "${CMAKE_CURRENT_SOURCE_DIR}/tinyexr/deps/miniz"
            PARENT_SCOPE)
        message(STATUS "Found RTX NTC target!")
    endif()
endif()
