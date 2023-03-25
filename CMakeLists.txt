
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

include_directories("/home/wanhsuan/OLSQ-dev/include/pybind11/include")

target_include_directories(olsq PUBLIC src include)

find_package(PythonInterp REQUIRED)
include_directories(${PYTHON_INCLUDE_DIRS})
# find_package(Python 3.6 COMPONENTS Interpreter Development REQUIRED)
# find_package(pybind11 CONFIG REQUIRED)
# find_package(pybind11 REQUIRED)
add_subdirectory("include/pybind11")
set(PYBIND11_CPP_STANDARD -std=c++14)

find_package(Z3 REQUIRED CONFIG)
set(Z3_DIR "/home/wanhsuan/z3/build")
set(PYBIND11_DIR "/home/wanhsuan/OLSQ-dev/include/pybind11")
find_package(Python REQUIRED COMPONENTS Development)

target_link_libraries(olsq PUBLIC z3::libz3 Python::Python)

add_executable(main app/main.cpp)
target_link_libraries(main olsq)

# pybind11_add_module(olsqPy ${SRC} ${API})
pybind11_add_module(olsqPy ${SRC} ${API})
target_include_directories(olsqPy PUBLIC src)
target_link_libraries(olsqPy PUBLIC z3::libz3)
