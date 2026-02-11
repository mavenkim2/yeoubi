find_package(pxr CONFIG REQUIRED HINTS "${USD_ROOT}")

if(pxr_FOUND)
    message(STATUS "Found OpenUSD: ${PXR_CMAKE_DIR}")
    set (OPENUSD_FOUND ON)
    set (OPENUSD_INCLUDE_DIRS ${PXR_INCLUDE_DIRS})
    set (OPENUSD_LIBRARIES "")

    set(TEMP_PXR_LIBRARY_DIR ${PXR_CMAKE_DIR}/lib)
    set(OPENUSD_LIBRARY_DIR ${PXR_CMAKE_DIR}/lib)

    foreach (LIB ${PXR_LIBRARIES})
        find_library(TEMP_LIB 
            NAMES "usd_${LIB}" "${LIB}" 
            PATHS ${TEMP_PXR_LIBRARY_DIR} 
            NO_DEFAULT_PATH
        )
        if (TEMP_LIB)
            list (APPEND OPENUSD_LIBRARIES ${TEMP_LIB})
        endif()
        unset (TEMP_LIB CACHE)
    endforeach()

    # Find TBB and OpenSubdiv libraries
    find_library(TEMP_TBB_LIBRARY_DEBUG_PXR NAMES tbb_debug tbb PATHS ${TEMP_PXR_LIBRARY_DIR} NO_DEFAULT_PATH)
    find_library(TEMP_TBB_LIBRARY_RELEASE_PXR NAMES tbb PATHS ${TEMP_PXR_LIBRARY_DIR} NO_DEFAULT_PATH)
    find_library(TEMP_OPENSUBDIV_LIBRARY_DEBUG_PXR NAMES osdCPU_d osdCPU PATHS ${TEMP_PXR_LIBRARY_DIR} NO_DEFAULT_PATH)
    find_library(TEMP_OPENSUBDIV_LIBRARY_RELEASE_PXR NAMES osdCPU PATHS ${TEMP_PXR_LIBRARY_DIR} NO_DEFAULT_PATH)

    if(TEMP_TBB_LIBRARY_RELEASE_PXR)
        # set(TBB_INCLUDE_DIRS ${PXR_INCLUDE_DIRS})
        list (APPEND OPENUSD_LIBRARIES 
                "$<$<CONFIG:Debug>:${TEMP_TBB_LIBRARY_DEBUG_PXR}>"
                "$<$<NOT:$<CONFIG:Debug>>:${TEMP_TBB_LIBRARY_RELEASE_PXR}>"
        )
        # set(TBB_LIBRARIES
        #     optimized ${_tbb_library_release_pxr}
        #     debug ${_tbb_library_debug_pxr}
        # )
        set(USD_OVERRIDE_TBB ON)
    endif()
    if (TEMP_OPENSUBDIV_LIBRARY_RELEASE_PXR)
        list (APPEND OPENUSD_LIBRARIES 
            "$<$<CONFIG:Debug>:${TEMP_OPENSUBDIV_LIBRARY_DEBUG_PXR}>"
            "$<$<NOT:$<CONFIG:Debug>>:${TEMP_OPENSUBDIV_LIBRARY_RELEASE_PXR}>"
        )
    endif()

    unset(TEMP_TBB_LIBRARY_DEBUG_PXR CACHE)
    unset(TEMP_TBB_LIBRARY_RELEASE_PXR CACHE)
    unset(TEMP_OPENSUBDIV_LIBRARY_DEBUG_PXR CACHE)
    unset(TEMP_OPENSUBDIV_LIBRARY_RELEASE_PXR CACHE)
    unset(TEMP_PXR_LIBRARY_DIR CACHE)
endif()

