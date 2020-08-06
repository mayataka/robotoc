cmake_minimum_required(VERSION 3.5)
project(iiwa14 CXX)

set(CMAKE_CXX_STANDARD 14)

find_package(idocp REQUIRED)

add_executable(
    ocp_test
    ocp_test.cpp
)
target_link_libraries(
    ocp_test
    PRIVATE
    idocp::idocp
)
target_include_directories(
    ocp_test
    PRIVATE
    ${IDOCP_INCLUDE_DIR}
)

# add_library(
#     mpc
#     SHARED 
#     ${IDOCP_SRC_DIR}/ocp/mpc.cpp
# )
# target_include_directories(
#     mpc
#     PRIVATE
#     ${EIGEN3_INCLUDE_DIR}
#     ${IDOCP_INCLUDE_DIR}
#     ${Boost_INCLUDE_DIRS}
#     ${PINOCCHIO_INCLUDEDIR}
# )
# target_link_libraries(
#     mpc
#     PRIVATE
#     idocp
#     robot
#     constraints
# )

# set(IDOCP_SIMULATOR_INCLUDE_DIR ${PROJECT_SOURCE_DIR}/../../simulator/include)
# set(IDOCP_SIMULATOR_SRC_DIR ${PROJECT_SOURCE_DIR}/../../simulator/src)

# add_library(
#     mpc_simulator
#     SHARED 
#     ${IDOCP_SIMULATOR_SRC_DIR}/simulator.cpp
#     ${IDOCP_SIMULATOR_SRC_DIR}/runge_kutta.cpp
#     ${IDOCP_SIMULATOR_SRC_DIR}/simulation_data_saver.cpp
# )
# target_include_directories(
#     mpc_simulator
#     PRIVATE
#     ${EIGEN3_INCLUDE_DIR}
#     ${IDOCP_INCLUDE_DIR}
#     ${Boost_INCLUDE_DIRS}
#     ${PINOCCHIO_INCLUDEDIR}
#     ${IDOCP_SIMULATOR_INCLUDE_DIR}
# )
# target_link_libraries( 
#     mpc_simulator
#     PRIVATE
#     mpc
#     idocp
#     robot
#     constraints
# )


# add_executable(
#     sim_mpc
#     sim_mpc.cpp
# )
# target_link_libraries(
#     sim_mpc
#     PRIVATE
#     mpc_simulator
#     mpc
#     idocp
#     robot
#     constraints
# )
# target_include_directories(
#     sim_mpc
#     PRIVATE
#     ${EIGEN3_INCLUDE_DIR}
#     ${IDOCP_INCLUDE_DIR}
#     ${Boost_INCLUDE_DIRS}
#     ${PINOCCHIO_INCLUDEDIR}
#     ${IDOCP_SIMULATOR_INCLUDE_DIR}
# )
