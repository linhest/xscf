
cmake_minimum_required (VERSION 2.6)
project (xscf)
find_package( BLAS REQUIRED )
find_package( LAPACK REQUIRED )
find_library ( LIBCINT_LIBRARY cint )
include_directories(${BLAS_INCLUDE_DIR})
set(CMAKE_BINARY_DIR ${CMAKE_SOURCE_DIR}/bin)
set(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR})
set(LIBRARY_OUTPUT_PATH ${CMAKE_BINARY_DIR})

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

#set(CMAKE_BUILD_TYPE Debug)
set(CMAKE_CXX_FLAGS_RELEASE "-O3")
set(CMAKE_CXX_FLAGS "-Wall -Wextra")
set(CMAKE_CXX_FLAGS_DEBUG "-g")

include_directories(include)
file(GLOB SOURCES "src/*.cc")
add_executable(xscf ${SOURCES})
target_link_libraries (xscf m blas lapack ${BLAS_LIBRARIES} ${LIBCINT_LIBRARY})

enable_testing()

add_test(
  NAME He
  COMMAND sh -c "cd ${CMAKE_SOURCE_DIR}; bin/xscf test_He/input.txt > tmp; diff tmp test_He/output.txt  | grep 'Energy=' || echo ok"
)
set_tests_properties(He PROPERTIES
                       PASS_REGULAR_EXPRESSION "ok")
add_test(
  NAME Ne
  COMMAND sh -c "cd ${CMAKE_SOURCE_DIR}; bin/xscf test_Ne/input.txt > tmp; diff tmp test_Ne/output.txt  | grep 'Energy=' || echo ok"
)
set_tests_properties(Ne PROPERTIES
                       PASS_REGULAR_EXPRESSION "ok")
add_test(
  NAME 2Ne
  COMMAND sh -c "cd ${CMAKE_SOURCE_DIR}; bin/xscf test_2Ne/input.txt > tmp; diff tmp test_2Ne/output.txt  | grep 'Energy=' || echo ok"
)
set_tests_properties(2Ne PROPERTIES
                       PASS_REGULAR_EXPRESSION "ok")
add_test(
  NAME H2O
  COMMAND sh -c "cd ${CMAKE_SOURCE_DIR}; bin/xscf test_H2O/input.txt > tmp; diff tmp test_H2O/output.txt  | grep 'Energy=' || echo ok"
)
set_tests_properties(H2O PROPERTIES
                       PASS_REGULAR_EXPRESSION "ok")
