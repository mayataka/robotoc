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
    idocp_robot
)

# build idocp constraints
add_library(
    idocp_constraints
    ${IDOCP_SRC_DIR}/constraints/pdipm_func.cpp
    ${IDOCP_SRC_DIR}/constraints/joint_space_constraints/joint_space_constraints.cpp
    ${IDOCP_SRC_DIR}/constraints/joint_space_constraints/joint_variables_upper_limits.cpp
    ${IDOCP_SRC_DIR}/constraints/joint_space_constraints/joint_variables_lower_limits.cpp
)
target_include_directories(
    idocp_constraints
    PUBLIC
    ${IDOCP_INCLUDE_DIR}
    ${Boost_INCLUDE_DIR}
    ${PINOCCHIO_INCLUDE_DIRS}
)
target_link_libraries(
    idocp_constraints
    PUBLIC
    idocp_robot
)

# build idocp manipulator
add_library(
    idocp_manipulator
    ${IDOCP_SRC_DIR}/manipulator/cost_function.cpp
    ${IDOCP_SRC_DIR}/manipulator/constraints.cpp
)
target_include_directories(
    idocp_manipulator
    PUBLIC
    ${IDOCP_INCLUDE_DIR}
    ${Boost_INCLUDE_DIR}
    ${PINOCCHIO_INCLUDE_DIRS}
)
target_link_libraries(
    idocp_manipulator
    PUBLIC
    idocp_cost
    idocp_robot
    idocp_constraints
)

# build idocp quadruped
add_library(
    idocp_quadruped
    ${IDOCP_SRC_DIR}/quadruped/cost_function.cpp
    ${IDOCP_SRC_DIR}/quadruped/constraints.cpp
)
target_include_directories(
    idocp_quadruped
    PUBLIC
    ${IDOCP_INCLUDE_DIR}
    ${Boost_INCLUDE_DIR}
    ${PINOCCHIO_INCLUDE_DIRS}
)
target_link_libraries(
    idocp_quadruped
    PUBLIC
    idocp_cost
    idocp_robot
    idocp_constraints
)

# build idocp ocp
add_library(
    idocp_ocp
    ${IDOCP_SRC_DIR}/ocp/riccati_matrix_factorizer.cpp
    ${IDOCP_SRC_DIR}/ocp/riccati_matrix_inverter.cpp
    ${IDOCP_SRC_DIR}/ocp/split_terminal_ocp.cpp
    # ${IDOCP_SRC_DIR}/ocp/split_ocp.cpp
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
    idocp_cost
    idocp_robot
    idocp_constraints
    idocp_manipulator
    idocp_quadruped
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