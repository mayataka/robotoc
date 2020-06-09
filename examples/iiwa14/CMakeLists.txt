cmake_minimum_required(VERSION 3.5)
project(iiwa14 CXX)

set(CMAKE_CXX_STANDARD 11)


set(IDOCP_INCLUDE_DIR ${PROJECT_SOURCE_DIR}/../../include/IDOCP)
set(IDOCP_SRC_DIR ${PROJECT_SOURCE_DIR}/../../src)
set(IIWA_DIR ${PROJECT_SOURCE_DIR})

find_package(Eigen3 REQUIRED)
find_package(Boost COMPONENTS system filesystem REQUIRED )
find_package(PkgConfig)
pkg_check_modules(PINOCCHIO REQUIRED pinocchio)
add_definitions(
    -DPINOCCHIO_URDFDOM_TYPEDEF_SHARED_PTR
    -DPINOCCHIO_URDFDOM_USE_STD_SHARED_PTR
    -DPINOCCHIO_WITH_URDFDOM
)
link_directories(${PINOCCHIO_LIBDIR})

find_package(OpenMP REQUIRED)
if(OpenMP_FOUND)
    set(CMAKE_C_FLAGS  "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

add_library(
    robot
    SHARED 
    ${IDOCP_SRC_DIR}/robot/robot.cpp
    ${IDOCP_SRC_DIR}/robot/passive_joints.cpp
    ${IDOCP_SRC_DIR}/robot/point_contact.cpp
)
target_include_directories(
    robot
    PRIVATE
    ${EIGEN3_INCLUDE_DIR}
    ${IDOCP_INCLUDE_DIR}
    ${Boost_INCLUDE_DIRS}
    ${PINOCCHIO_INCLUDEDIR}
)
target_link_libraries(
    robot
    PRIVATE
    pthread
    ${Boost_LIBRARIES}
    urdfdom_model
    ${PINOCCHIO_LIBRARIES}
)

add_library(
    invdynocp
    SHARED 
    ${IDOCP_SRC_DIR}/ocp/OCP.cpp
    ${IDOCP_SRC_DIR}/ocp/split_OCP.cpp
    ${IDOCP_SRC_DIR}/cost/configuration_space_cost.cpp
    ${IDOCP_SRC_DIR}/constraints/configuration_space_constraints.cpp
)
target_include_directories(
    invdynocp
    PRIVATE
    ${EIGEN3_INCLUDE_DIR}
    ${IDOCP_INCLUDE_DIR}
    ${Boost_INCLUDE_DIRS}
    ${PINOCCHIO_INCLUDEDIR}
)
target_link_libraries(
    invdynocp
    PRIVATE
    robot
)

add_executable(
    a.out 
    main.cpp
    cost_function.cpp
    constraints.cpp
)
target_link_libraries(
    a.out
    PRIVATE
    robot
    invdynocp
)
target_include_directories(
    a.out 
    PRIVATE
    ${EIGEN3_INCLUDE_DIR}
    ${IDOCP_INCLUDE_DIR}
    ${Boost_INCLUDE_DIRS}
    ${PINOCCHIO_INCLUDEDIR}
)
target_compile_options(
    a.out
    PRIVATE
    -O3
)