cmake_minimum_required(VERSION 3.5)
project(iiwa14 CXX)

set(CMAKE_CXX_STANDARD 11)


set(IDOCP_INCLUDE_DIR ${PROJECT_SOURCE_DIR}/../../include/idocp)
set(IDOCP_SRC_DIR ${PROJECT_SOURCE_DIR}/../../src)

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
    ${IDOCP_SRC_DIR}/robot/floating_base.cpp
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
    constraints
    SHARED 
    ${IDOCP_SRC_DIR}/constraints/joint_space_constraints/joint_space_constraints.cpp
    ${IDOCP_SRC_DIR}/constraints/joint_space_constraints/joint_position_upper_limits.cpp
    ${IDOCP_SRC_DIR}/constraints/joint_space_constraints/joint_position_lower_limits.cpp
    ${IDOCP_SRC_DIR}/constraints/joint_space_constraints/joint_velocity_upper_limits.cpp
    ${IDOCP_SRC_DIR}/constraints/joint_space_constraints/joint_velocity_lower_limits.cpp
    ${IDOCP_SRC_DIR}/constraints/joint_space_constraints/joint_torque_upper_limits.cpp
    ${IDOCP_SRC_DIR}/constraints/joint_space_constraints/joint_torque_lower_limits.cpp
    ${IDOCP_SRC_DIR}/constraints/pdipm_func.cpp
)
target_include_directories(
    constraints
    PRIVATE
    ${EIGEN3_INCLUDE_DIR}
    ${IDOCP_INCLUDE_DIR}
    ${Boost_INCLUDE_DIRS}
    ${PINOCCHIO_INCLUDEDIR}
)
target_link_libraries(
    constraints
    PRIVATE
    robot
)


add_library(
    idocp
    SHARED 
    ${IDOCP_SRC_DIR}/ocp/ocp.cpp
    ${IDOCP_SRC_DIR}/ocp/split_ocp.cpp
    ${IDOCP_SRC_DIR}/ocp/line_search_filter.cpp
    ${IDOCP_SRC_DIR}/cost/joint_space_cost.cpp
    ${IDOCP_SRC_DIR}/cost/contact_cost.cpp
    ${IDOCP_SRC_DIR}/manipulator/cost_function.cpp
    ${IDOCP_SRC_DIR}/manipulator/constraints.cpp
)
target_include_directories(
    idocp
    PRIVATE
    ${EIGEN3_INCLUDE_DIR}
    ${IDOCP_INCLUDE_DIR}
    ${Boost_INCLUDE_DIRS}
    ${PINOCCHIO_INCLUDEDIR}
)
target_link_libraries(
    idocp
    PRIVATE
    robot
    constraints
)


add_executable(
    ocp_test
    ocp_test.cpp
)
target_link_libraries(
    ocp_test
    PRIVATE
    robot
    idocp
)
target_include_directories(
    ocp_test
    PRIVATE
    ${EIGEN3_INCLUDE_DIR}
    ${IDOCP_INCLUDE_DIR}
    ${Boost_INCLUDE_DIRS}
    ${PINOCCHIO_INCLUDEDIR}
)

add_library(
    mpc
    SHARED 
    ${IDOCP_SRC_DIR}/ocp/mpc.cpp
)
target_include_directories(
    mpc
    PRIVATE
    ${EIGEN3_INCLUDE_DIR}
    ${IDOCP_INCLUDE_DIR}
    ${Boost_INCLUDE_DIRS}
    ${PINOCCHIO_INCLUDEDIR}
)
target_link_libraries(
    mpc
    PRIVATE
    idocp
    robot
    constraints
)

set(IDOCP_SIMULATOR_INCLUDE_DIR ${PROJECT_SOURCE_DIR}/../../simulator/include)
set(IDOCP_SIMULATOR_SRC_DIR ${PROJECT_SOURCE_DIR}/../../simulator/src)

add_library(
    mpc_simulator
    SHARED 
    ${IDOCP_SIMULATOR_SRC_DIR}/simulator.cpp
    ${IDOCP_SIMULATOR_SRC_DIR}/runge_kutta.cpp
    ${IDOCP_SIMULATOR_SRC_DIR}/simulation_data_saver.cpp
)
target_include_directories(
    mpc_simulator
    PRIVATE
    ${EIGEN3_INCLUDE_DIR}
    ${IDOCP_INCLUDE_DIR}
    ${Boost_INCLUDE_DIRS}
    ${PINOCCHIO_INCLUDEDIR}
    ${IDOCP_SIMULATOR_INCLUDE_DIR}
)
target_link_libraries( 
    mpc_simulator
    PRIVATE
    mpc
    idocp
    robot
    constraints
)


add_executable(
    sim_mpc
    sim_mpc.cpp
)
target_link_libraries(
    sim_mpc
    PRIVATE
    mpc_simulator
    mpc
    idocp
    robot
    constraints
)
target_include_directories(
    sim_mpc
    PRIVATE
    ${EIGEN3_INCLUDE_DIR}
    ${IDOCP_INCLUDE_DIR}
    ${Boost_INCLUDE_DIRS}
    ${PINOCCHIO_INCLUDEDIR}
    ${IDOCP_SIMULATOR_INCLUDE_DIR}
)
