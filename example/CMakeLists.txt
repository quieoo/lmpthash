cmake_minimum_required(VERSION 3.10)

project(HiDPU_Simulator)

set(CMAKE_C_STANDARD 99)

set(CMAKE_CXX_STANDARD 11)

set(SOURCE_FILES HiDPU_Simulator.c)

add_executable(hidpu ${SOURCE_FILES})

find_library(CLMPTHASH_LIB clmpthash PATHS ${CMAKE_SOURCE_DIR}/../build)

if (CLMPTHASH_LIB)
    message(STATUS "Found clmpthash library: ${CLMPTHASH_LIB}")
    target_link_libraries(hidpu ${CLMPTHASH_LIB} stdc++ m pthread)
else()
    message(FATAL_ERROR "clmpthash library not found")
endif()


target_link_libraries(hidpu ${CLMPTHASH_LIB} stdc++ m pthread )
