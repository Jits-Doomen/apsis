cmake_minimum_required(VERSION 3.20)
project(Apsis VERSION 0.2.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

add_library(apsis_core
    src/interpreter.cpp
    packages/formulas.cpp
    packages/ephemerides.cpp
    packages/propagator.cpp
    packages/trajectory.cpp
    src/numerics.cpp
    src/linear_algebra.cpp
    src/orbit.cpp
    src/units.cpp
    src/data_io.cpp
)

target_include_directories(apsis_core PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/include")

if(MSVC)
  target_compile_options(apsis_core PRIVATE /W4 /permissive-)
else()
  target_compile_options(apsis_core PRIVATE -Wall -Wextra -Wpedantic)
endif()

add_executable(apsis main.cpp)
target_link_libraries(apsis PRIVATE apsis_core)

include(CTest)
if(BUILD_TESTING)
  add_executable(apsis_science_tests tests/test_science.cpp)
  target_link_libraries(apsis_science_tests PRIVATE apsis_core)
  add_test(NAME scientific_core COMMAND apsis_science_tests)
endif()

install(TARGETS apsis RUNTIME DESTINATION bin)