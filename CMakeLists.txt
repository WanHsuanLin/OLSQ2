cmake_minimum_required(VERSION 3.15)
project(OLSQ VERSION 0.1.0 LANGUAGES CXX)

# global config
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# dependencies
find_package(Python REQUIRED COMPONENTS Interpreter Development)
find_package(Z3 REQUIRED)
add_subdirectory(pybind11)

# olsq library
add_library(olsq
  src/cir/circuit.cpp
  src/device/device.cpp
  src/misc/timeUsage.cpp
  src/olsq/olsq.cpp
)

set_target_properties(olsq PROPERTIES POSITION_INDEPENDENT_CODE TRUE)

target_include_directories(olsq
  PUBLIC
    "$<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/src>"
    "$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>"
)

target_link_libraries(olsq PUBLIC z3::libz3)

# python bindings
pybind11_add_module(olsqPy src/api/apiPy.cpp)
target_link_libraries(olsqPy PUBLIC olsq)

# main program
add_executable(main app/main.cpp)
target_link_libraries(main olsq)
