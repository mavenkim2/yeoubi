# NOTE: this was build by running the provided build script in OpenUSD with --no-python,
# then setting USD_ROOT to the install location.
set(USD_ROOT "" CACHE PATH "Path to the OpenUSD installation directory")
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

set(OPENCOLORIO_ROOT "${CMAKE_CURRENT_SOURCE_DIR}/OpenColorIO")
if (EXISTS "${OPENCOLORIO_ROOT}/CMakeLists.txt")
    set(BUILD_SHARED_LIBS ON CACHE BOOL "" FORCE)
    set(OCIO_BUILD_APPS OFF CACHE BOOL "" FORCE)
    set(OCIO_BUILD_TESTS OFF CACHE BOOL "" FORCE)
    set(OCIO_BUILD_GPU_TESTS OFF CACHE BOOL "" FORCE)
    set(OCIO_BUILD_DOCS OFF CACHE BOOL "" FORCE)
    set(OCIO_BUILD_PYTHON OFF CACHE BOOL "" FORCE)
    set(OCIO_BUILD_JAVA OFF CACHE BOOL "" FORCE)
    set(OCIO_INSTALL_EXT_PACKAGES MISSING CACHE STRING "" FORCE)
    add_subdirectory(OpenColorIO EXCLUDE_FROM_ALL)

    if (TARGET OpenColorIO)
        set(OPENCOLORIO_FOUND ON PARENT_SCOPE)
        set(OPENCOLORIO_INCLUDE_DIRS "${OPENCOLORIO_ROOT}/include" PARENT_SCOPE)
        set(OPENCOLORIO_LIBRARIES OpenColorIO PARENT_SCOPE)
        message(STATUS "Found OpenColorIO target!")
    endif()
endif()

add_subdirectory(meshoptimizer)

if (TARGET meshoptimizer)
    set (MESHOPTIMIZER_FOUND ON PARENT_SCOPE)
    set(MESHOPTIMIZER_INCLUDE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/meshoptimizer/src" PARENT_SCOPE)
    set(MESHOPTIMIZER_LIBRARY meshoptimizer PARENT_SCOPE)
    message(STATUS "Found meshoptimizer target!")
endif()

