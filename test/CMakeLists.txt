cmake_minimum_required(VERSION 3.1)
project(idocp_test CXX)
set(CMAKE_CXX_STANDARD 11)

# enable coverage
set(CMAKE_C_FLAGS  "${CMAKE_C_FLAGS} --coverage")
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} --coverage")

# enable google test
enable_testing()
include(cmake/gtest.cmake)

# find Eigen
find_package(Eigen3 REQUIRED)
include_directories(${EIGEN3_INCLUDE_DIR})

# find Boost
find_package(Boost COMPONENTS system filesystem REQUIRED)

# find pinocchio
find_package(PkgConfig)
pkg_check_modules(PINOCCHIO REQUIRED pinocchio)
add_definitions(
    -DPINOCCHIO_URDFDOM_TYPEDEF_SHARED_PTR
    -DPINOCCHIO_URDFDOM_USE_STD_SHARED_PTR
    -DPINOCCHIO_WITH_URDFDOM
)
link_directories(${PINOCCHIO_LIBDIR})

# set IDOCP directories
set(IDOCP_SRC_DIR ${CMAKE_CURRENT_SOURCE_DIR}/../src)
set(IDOCP_INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/../include/idocp)

# coverage options
add_library(coverage_config INTERFACE)
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
  target_compile_options(coverage_config INTERFACE
    -O0        
    -g         
    --coverage 
  )
  if(CMAKE_VERSION VERSION_GREATER_EQUAL 3.13)
    target_link_options(coverage_config INTERFACE --coverage)
  else()
    target_link_libraries(coverage_config INTERFACE --coverage)
  endif()
endif()

add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/robot)
# add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/constraints)