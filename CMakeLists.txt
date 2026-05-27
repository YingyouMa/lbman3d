cmake_minimum_required(VERSION 3.23)

project(LBM)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_FLAGS "-fopenmp -g -Ofast")

# Detect whether this toolchain's standard library supports <format> and <print>.
# GCC 13's libstdc++ ships incomplete C++23 headers; GCC 15 is fine.
include(CheckSourceCompiles)
set(CMAKE_REQUIRED_FLAGS "-std=c++23")
check_source_compiles(CXX
    "#include <format>\n#include <print>\nint main(){std::println(\"{}\",42);}"
    HAVE_STD_FORMAT
)

include(FetchContent)

# Fetch Kokkos mdspan
FetchContent_Declare(
        mdspan
        GIT_REPOSITORY https://github.com/kokkos/mdspan.git
        GIT_TAG stable
)
FetchContent_MakeAvailable(mdspan)

# Fetch {fmt} only when the standard library lacks <format>/<print>
if(NOT HAVE_STD_FORMAT)
    FetchContent_Declare(
        fmt
        GIT_REPOSITORY https://github.com/fmtlib/fmt.git
        GIT_TAG 12.1.0
    )
    FetchContent_MakeAvailable(fmt)
endif()

find_package(HDF5 REQUIRED COMPONENTS C)

add_executable(main)
set_target_properties(main PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})

target_include_directories(main PUBLIC
        ${HDF5_INCLUDE_DIRS}
)

if(DEFINED LBMAN3D_CASE_DIR)
    target_include_directories(main BEFORE PUBLIC
        ${LBMAN3D_CASE_DIR}
    )
endif()

target_link_libraries(main PRIVATE
    ${HDF5_C_LIBRARIES})

target_include_directories(main PUBLIC
        ${mdspan_SOURCE_DIR}
)

target_link_libraries(main PRIVATE mdspan)

# Wire up std::format detection result into the build
if(HAVE_STD_FORMAT)
    target_compile_definitions(main PRIVATE HAVE_STD_FORMAT)
else()
    target_link_libraries(main PRIVATE fmt::fmt)
endif()

add_subdirectory(src)
