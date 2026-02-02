# NOTE: this was build by running the provided build script in OpenUSD with --no-python,
# then setting USD_ROOT to the install location.
set(USD_ROOT "" CACHE PATH "Path to the OpenUSD installation directory")

set (OLD_CMAKE_PREFIX_PATH ${CMAKE_PREFIX_PATH})
list(APPEND CMAKE_PREFIX_PATH "${USD_ROOT}/lib/cmake")
find_package(pxr HINTS "${USD_ROOT}")

if(pxr_FOUND)
    message(STATUS "Found OpenUSD version: ${pxr_VERSION}")
    set (OPENUSD_LIBRARIES usd usdGeom PARENT_SCOPE)
    set (OPENUSD_INCLUDE_DIRS ${PXR_INCLUDE_DIRS} PARENT_SCOPE)
endif()

set (CMAKE_PREFIX_PATH ${OLD_CMAKE_PREFIX_PATH})
