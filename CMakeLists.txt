cmake_minimum_required(VERSION 3.15)
project(lx-mie VERSION 1.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

add_compile_options(-Wall -pedantic)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR})

add_subdirectory(src/mie)
add_subdirectory(src/main)

option(BUILD_PYTHON_BINDINGS "Build Python bindings via PyBind11" OFF)
if(BUILD_PYTHON_BINDINGS OR SKBUILD)
  add_subdirectory(src/python)
endif()
