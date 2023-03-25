
# find_package(PythonInterp REQUIRED)
# add_subdirectory(include/pybind11)

# pybind11_add_module(olsq ${SOURCE_FILES})

cmake_minimum_required(VERSION 3.1)
project(OLSQ
VERSION 0.1.0
LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

add_library(olsq
  src/cir/circuit.cpp
  src/device/device.cpp
  src/misc/timeUsage.cpp
  src/olsq/olsq.cpp
  src/api/apiPy.cpp
)

set(SRC
  src/cir/circuit.cpp
  src/device/device.cpp
  src/misc/timeUsage.cpp
  src/olsq/olsq.cpp
  src/api/apiPy.cpp
)

set(API
   src/api/apiPy.cpp
 )


target_include_directories(olsq PUBLIC src include)

find_package(PythonInterp REQUIRED)
include_directories(${PYTHON_INCLUDE_DIRS})

# Find pybind11
if (PYBIND11_DIR)
  set(PYBIND11_ROOT_DIR ${PYBIND11_DIR})
else()
  set(PYBIND11_ROOT_DIR $ENV{PYBIND11_DIR})
endif()
message(STATUS "Find Pybind11: ${PYBIND11_ROOT_DIR}")
unset(PYBIND11_DIR)
# include_directories(${PYBIND11_ROOT_DIR})
# add_subdirectory(${PYBIND11_ROOT_DIR} "./pybind11")
include_directories("${PYBIND11_ROOT_DIR}/include")
add_subdirectory("include/pybind11")
set(PYBIND11_CPP_STANDARD -std=c++14)

find_package(Z3 REQUIRED CONFIG)

target_link_libraries(olsq PUBLIC z3::libz3 Python::Python)

add_executable(main app/main.cpp)
target_link_libraries(main olsq)

# pybind11_add_module(olsqPy ${SRC} ${API})
pybind11_add_module(olsqPy ${SRC} ${API})
target_include_directories(olsqPy PUBLIC src)
target_link_libraries(olsqPy PUBLIC z3::libz3)

# cmake .. -DPYBIND11_DIR=<Path_to_pybind11> -DZ3_DIR=<Path_to_z3_build>