cmake_minimum_required(VERSION 3.22.1)
if(NOT DEFINED CMAKE_TOOLCHAIN_FILE AND EXISTS "${CMAKE_BINARY_DIR}/conan_toolchain.cmake")
    set(CMAKE_TOOLCHAIN_FILE "${CMAKE_BINARY_DIR}/conan_toolchain.cmake")
endif()
project(JSSP_solver CXX)

set(CMAKE_CXX_STANDARD 17)

enable_testing()

set(ortools_DIR )

find_package(ortools CONFIG REQUIRED)
find_package(ZLIB CONFIG REQUIRED)

#include_directories(test/include)

add_executable(JSSP_solver main.cpp
        GeneticAlgorithm.cpp
        GeneticAlgorithm.h
        InstanceHandler.h
        Structs.h
        FunctionFactory.h
        Algorithm.h
        OrtoolsAlgorithm.h
        OrtoolsAlgorithm.cpp
        AlgorithmManager.h
        AlgorithmManager.cpp
        include/jssp.cpp
        include/jssp.h
        include/tabu_list.h
        include/ts.cpp
        include/ts.h
        include/json.hpp
        Logger.h
)

target_link_libraries(JSSP_solver PRIVATE ortools::ortools)

add_test(NAME run_test_ortools COMMAND test_ortools)