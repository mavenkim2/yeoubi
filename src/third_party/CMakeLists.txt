# NOTE: this was build by running the provided build script in OpenUSD with --no-python,
# then setting USD_ROOT to the install location.
set(USD_ROOT "" CACHE PATH "Path to the OpenUSD installation directory")

if (USD_ROOT)
    find_package(OpenUSD)
    if (OPENUSD_FOUND)
        set (WITH_USD ON PARENT_SCOPE) 
        set (OPENUSD_INCLUDE_DIRS ${OPENUSD_INCLUDE_DIRS} PARENT_SCOPE)
        set (OPENUSD_LIBRARIES ${OPENUSD_LIBRARIES} PARENT_SCOPE)
    endif()
endif()
