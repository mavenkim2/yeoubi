include(FindPackageHandleStandardArgs)

set(LIBNTC_ROOT "" CACHE PATH "Path to a prebuilt LibNTC installation")
set(LIBNTC_INCLUDE_DIR "" CACHE PATH "Path to the LibNTC include directory")
set(LIBNTC_LIBRARY "" CACHE FILEPATH "Path to the prebuilt libntc library")

set(_RTXNTC_ROOT_HINTS
    ${LIBNTC_ROOT}
    $ENV{LIBNTC_ROOT})

if (LIBNTC_INCLUDE_DIR)
    find_path(RTXNTC_INCLUDE_DIR
        NAMES libntc/ntc.h
        HINTS "${LIBNTC_INCLUDE_DIR}"
        NO_DEFAULT_PATH)
else()
    find_path(RTXNTC_INCLUDE_DIR
        NAMES libntc/ntc.h
        HINTS ${_RTXNTC_ROOT_HINTS}
        PATH_SUFFIXES include)
endif()

if (LIBNTC_LIBRARY)
    if (EXISTS "${LIBNTC_LIBRARY}" AND NOT IS_DIRECTORY "${LIBNTC_LIBRARY}")
        set(RTXNTC_LIBRARY "${LIBNTC_LIBRARY}")
    else()
        set(RTXNTC_LIBRARY "RTXNTC_LIBRARY-NOTFOUND")
    endif()
else()
    find_library(RTXNTC_LIBRARY
        NAMES libntc ntc
        HINTS ${_RTXNTC_ROOT_HINTS}
        PATH_SUFFIXES lib lib64 bin)
endif()

find_package_handle_standard_args(RTXNTC
    REQUIRED_VARS RTXNTC_INCLUDE_DIR RTXNTC_LIBRARY)

if (RTXNTC_FOUND)
    set(RTXNTC_INCLUDE_DIRS "${RTXNTC_INCLUDE_DIR}")
    set(RTXNTC_LIBRARIES libntc)

    if (NOT TARGET libntc)
        add_library(libntc UNKNOWN IMPORTED GLOBAL)
        set_target_properties(libntc PROPERTIES
            IMPORTED_LOCATION "${RTXNTC_LIBRARY}"
            INTERFACE_INCLUDE_DIRECTORIES "${RTXNTC_INCLUDE_DIR}")
    endif()
endif()

mark_as_advanced(RTXNTC_INCLUDE_DIR RTXNTC_LIBRARY)
