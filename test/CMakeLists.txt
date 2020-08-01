cmake_minimum_required(VERSION 3.1)
project(idocp_test CXX)
set(CMAKE_CXX_STANDARD 11)

# enable coverage
set(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} --coverage")
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

# set idocp directories
set(IDOCP_SRC_DIR ${CMAKE_CURRENT_SOURCE_DIR}/../src)
set(IDOCP_INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/../include/idocp)

# build idocp robot
add_library(
    idocp_robot
    ${IDOCP_SRC_DIR}/robot/robot.cpp
    ${IDOCP_SRC_DIR}/robot/floating_base.cpp
    ${IDOCP_SRC_DIR}/robot/point_contact.cpp
)
target_include_directories(
    idocp_robot
    PUBLIC
    ${IDOCP_INCLUDE_DIR}
    ${Boost_INCLUDE_DIR}
    ${PINOCCHIO_INCLUDE_DIRS}
)
target_link_libraries(
    idocp_robot
    PUBLIC
    pthread
    ${Boost_LIBRARIES}
    ${PINOCCHIO_LIBRARIES}
    urdfdom_model
)

# build idocp cost
add_library(
    idocp_cost
    ${IDOCP_SRC_DIR}/cost/joint_space_cost.cpp
    ${IDOCP_SRC_DIR}/cost/contact_cost.cpp
)
target_include_directories(
    idocp_cost
    PUBLIC
    ${IDOCP_INCLUDE_DIR}
    ${Boost_INCLUDE_DIR}
    ${PINOCCHIO_INCLUDE_DIRS}
)
target_link_libraries(
    idocp_cost
    PUBLIC
    robot
)

# build idocp ocp
add_library(
    idocp_ocp
    ${IDOCP_SRC_DIR}/ocp/riccati_matrix_factorizer.cpp
    ${IDOCP_SRC_DIR}/ocp/riccati_matrix_inverter.cpp
)
target_include_directories(
    idocp_ocp
    PUBLIC
    ${IDOCP_INCLUDE_DIR}
    ${Boost_INCLUDE_DIR}
    ${PINOCCHIO_INCLUDE_DIRS}
)
target_link_libraries(
    idocp_ocp
    PUBLIC
    cost
    robot
)

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
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/cost)
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/ocp)