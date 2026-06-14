cmake_minimum_required(VERSION 3.16)

project(Yeoubi LANGUAGES CXX) 

set(CMAKE_CXX_STANDARD 17 CACHE STRING "C++ standard for yeoubi and subprojects" FORCE)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/src/cmake")

set(YEOUBI_COMPILE_DEFINITIONS)
if (WIN32)
    list(APPEND YEOUBI_COMPILE_DEFINITIONS NOMINMAX)
endif()

add_executable(yeoubi
               ${CMAKE_CURRENT_SOURCE_DIR}/src/core/allocation.cpp
               ${CMAKE_CURRENT_SOURCE_DIR}/src/io/usd.cpp)

set(USD_ROOT "" CACHE PATH "Path to the OpenUSD installation directory")
if (USD_ROOT)
    find_package(OpenUSD)
    if (OPENUSD_FOUND)
        set (WITH_USD ON)
        set (OPENUSD_INCLUDE_DIRS ${OPENUSD_INCLUDE_DIRS})
        set (OPENUSD_LIBRARIES ${OPENUSD_LIBRARIES})
        set (OPENUSD_LIBRARY_DIR ${OPENUSD_LIBRARY_DIR})
    endif()
endif()

target_include_directories(yeoubi PRIVATE 
    ${CMAKE_CURRENT_SOURCE_DIR}/src
    ${OPENUSD_INCLUDE_DIRS}
)

target_link_libraries(yeoubi PRIVATE 
    ${OPENUSD_LIBRARIES}
)

target_compile_definitions(yeoubi PRIVATE
    ${YEOUBI_COMPILE_DEFINITIONS}
    $<$<CONFIG:Debug,RelWithDebInfo>:YBI_DEBUG>
)
