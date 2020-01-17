cmake_minimum_required(VERSION 3.15)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")

set(GCC_CLANG_WARNINGS
        -Wall -Wextra -Wshadow -Wnon-virtual-dtor -Wpedantic
        -Wold-style-cast -Wcast-align -Wunused -Woverloaded-virtual
        -Wconversion -Wsign-conversion -Wnull-dereference
        -Wdouble-promotion -Wformat=2)

set(GCC_EXTRA_WARNINGS
        -Wmisleading-indentation -Wduplicated-cond -Wduplicated-branches
        -Wlogical-op -Wuseless-cast)

project(exact_hubbard)

add_executable(exact_hubbard
        src/main.cpp
        src/always_false.hpp
        src/check_config.cpp
        src/io.hpp
        src/io.cpp
        src/linalg.hpp
        src/operator.hpp
        src/state.hpp
        src/state.cpp
        src/spectrum.cpp
        src/spectrum.hpp
        src/correlators.cpp
        src/correlators.hpp)

set_target_properties(exact_hubbard PROPERTIES
        CXX_STANDARD 17
        CXX_STANDARD_REQUIRED ON)

if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
    target_compile_options(exact_hubbard PUBLIC ${GCC_CLANG_WARNINGS})
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
    target_compile_options(exact_hubbard PUBLIC ${GCC_CLANG_WARNINGS} ${GCC_EXTRA_WARNINGS})
endif ()

if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang"
        OR "${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
    target_compile_options(exact_hubbard PUBLIC -march=native)
endif ()

target_link_libraries(exact_hubbard stdc++fs)

find_package(blaze REQUIRED)
target_include_directories(exact_hubbard SYSTEM PUBLIC ${blaze_INCLUDE_DIRS})
target_compile_options(exact_hubbard PUBLIC "${blaze_CXX_FLAGS}")
target_link_libraries(exact_hubbard ${blaze_LIBRARIES})
target_link_options(exact_hubbard PUBLIC "${blaze_LINKER_FLAGS}")
