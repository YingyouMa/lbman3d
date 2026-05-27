cmake_minimum_required(VERSION 3.23)
project(lbman3d_case)

if(EXISTS "${CMAKE_CURRENT_LIST_DIR}/../src/CMakeLists.txt")
    set(LBMAN3D_SOURCE_ROOT "${CMAKE_CURRENT_LIST_DIR}/..")
elseif(EXISTS "${CMAKE_CURRENT_LIST_DIR}/../../src/CMakeLists.txt")
    set(LBMAN3D_SOURCE_ROOT "${CMAKE_CURRENT_LIST_DIR}/../..")
else()
    message(FATAL_ERROR "Could not locate lbman3d source root relative to ${CMAKE_CURRENT_LIST_DIR}")
endif()

set(LBMAN3D_CASE_DIR "${CMAKE_CURRENT_LIST_DIR}" CACHE PATH "Case directory containing case_config_local.h" FORCE)
add_subdirectory("${LBMAN3D_SOURCE_ROOT}" lbman3d-build)
