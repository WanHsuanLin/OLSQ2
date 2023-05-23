cmake_minimum_required(VERSION 3.10)
project(OLSQ VERSION 0.1.0 LANGUAGES CXX)

# global config
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

add_subdirectory(pblib/include)

# dependencies
# find_package(Python REQUIRED COMPONENTS Interpreter Development)
find_package(PythonInterp REQUIRED)
find_package(PythonLibs REQUIRED)
include_directories(${PYTHON_INCLUDE_DIRS})
add_subdirectory(pybind11)

find_package(Bitwuzla REQUIRED)
find_package(PBlib REQUIRED)

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

target_link_libraries(olsq PUBLIC pblib/lib/libpb.a Bitwuzla::bitwuzla)

# python bindings
pybind11_add_module(olsqpyb src/api/apiPy.cpp)
target_link_libraries(olsqpyb PUBLIC olsq)

# python api
add_custom_target(pyolsq ALL
  COMMAND
    ${CMAKE_COMMAND} -E create_symlink "${PROJECT_SOURCE_DIR}/src/pyolsq/"  "${PROJECT_BINARY_DIR}/pyolsq"
)

if(CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME)
  include(CTest)
endif()

if(CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME AND BUILD_TESTING)
  add_executable(simple test/simple.cpp)
  target_link_libraries(simple olsq)
  add_test(NAME simple-test COMMAND simple)
endif()


