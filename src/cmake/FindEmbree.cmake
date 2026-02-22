include(FindPackageHandleStandardArgs)

set(_EMBREE_ROOT_HINTS
    ${EMBREE_ROOT}
    $ENV{EMBREE_ROOT}
    $ENV{EMBREE_LOCATION})

find_path(EMBREE_INCLUDE_DIR
    NAMES embree4/rtcore.h embree3/rtcore.h
    HINTS ${_EMBREE_ROOT_HINTS}
    PATH_SUFFIXES include)

find_library(EMBREE_LIBRARY
    NAMES embree4 embree3 embree
    HINTS ${_EMBREE_ROOT_HINTS}
    PATH_SUFFIXES lib lib64)

set(EMBREE_VERSION "")
set(EMBREE_VERSION_MAJOR "")
set(EMBREE_VERSION_MINOR "")
set(EMBREE_VERSION_PATCH "")

set(_EMBREE_RTCORE_CONFIG "")
if (EMBREE_INCLUDE_DIR)
    if (EXISTS "${EMBREE_INCLUDE_DIR}/embree4/rtcore_config.h")
        set(_EMBREE_RTCORE_CONFIG "${EMBREE_INCLUDE_DIR}/embree4/rtcore_config.h")
        set(EMBREE_VERSION_MAJOR 4)
    elseif (EXISTS "${EMBREE_INCLUDE_DIR}/embree3/rtcore_config.h")
        set(_EMBREE_RTCORE_CONFIG "${EMBREE_INCLUDE_DIR}/embree3/rtcore_config.h")
        set(EMBREE_VERSION_MAJOR 3)
    endif()
endif()

if (_EMBREE_RTCORE_CONFIG)
    file(STRINGS "${_EMBREE_RTCORE_CONFIG}" _EMBREE_MAJOR_LINE
         REGEX "^#define RTC_VERSION_MAJOR[ \t]+[0-9]+")
    file(STRINGS "${_EMBREE_RTCORE_CONFIG}" _EMBREE_MINOR_LINE
         REGEX "^#define RTC_VERSION_MINOR[ \t]+[0-9]+")
    file(STRINGS "${_EMBREE_RTCORE_CONFIG}" _EMBREE_PATCH_LINE
         REGEX "^#define RTC_VERSION_PATCH[ \t]+[0-9]+")

    if (_EMBREE_MAJOR_LINE)
        string(REGEX REPLACE ".*RTC_VERSION_MAJOR[ \t]+([0-9]+).*" "\\1"
               EMBREE_VERSION_MAJOR "${_EMBREE_MAJOR_LINE}")
    endif()
    if (_EMBREE_MINOR_LINE)
        string(REGEX REPLACE ".*RTC_VERSION_MINOR[ \t]+([0-9]+).*" "\\1"
               EMBREE_VERSION_MINOR "${_EMBREE_MINOR_LINE}")
    endif()
    if (_EMBREE_PATCH_LINE)
        string(REGEX REPLACE ".*RTC_VERSION_PATCH[ \t]+([0-9]+).*" "\\1"
               EMBREE_VERSION_PATCH "${_EMBREE_PATCH_LINE}")
    endif()
endif()

if (EMBREE_VERSION_MAJOR)
    if (EMBREE_VERSION_MINOR AND EMBREE_VERSION_PATCH)
        set(EMBREE_VERSION
            "${EMBREE_VERSION_MAJOR}.${EMBREE_VERSION_MINOR}.${EMBREE_VERSION_PATCH}")
    else()
        set(EMBREE_VERSION "${EMBREE_VERSION_MAJOR}")
    endif()
endif()

find_package_handle_standard_args(Embree
    REQUIRED_VARS EMBREE_INCLUDE_DIR EMBREE_LIBRARY
    VERSION_VAR EMBREE_VERSION)

if (EMBREE_FOUND)
    set(EMBREE_INCLUDE_DIRS "${EMBREE_INCLUDE_DIR}")
    set(EMBREE_LIBRARIES "${EMBREE_LIBRARY}")
endif()

mark_as_advanced(EMBREE_INCLUDE_DIR EMBREE_LIBRARY)

