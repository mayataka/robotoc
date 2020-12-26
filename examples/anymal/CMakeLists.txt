cmake_minimum_required(VERSION 3.5)
project(anymal CXX)

set(CMAKE_CXX_STANDARD 14)

find_package(idocp REQUIRED)
find_package(raisim CONFIG REQUIRED)
find_package(raisimOgre CONFIG REQUIRED)

add_executable(
    ocp_benchmark
    ocp_benchmark.cpp
)
target_link_libraries(
    ocp_benchmark
    PRIVATE
    idocp::idocp
    raisim::raisim
    raisim::raisimOgre
)
target_include_directories(
    ocp_benchmark
    PRIVATE
    ${IDOCP_INCLUDE_DIR}
)

add_executable(
    ocp_trotting
    ocp_trotting.cpp
)
target_link_libraries(
    ocp_trotting
    PRIVATE
    idocp::idocp
    raisim::raisim
    raisim::raisimOgre
)
target_include_directories(
    ocp_trotting
    PRIVATE
    ${IDOCP_INCLUDE_DIR}
)

option(SIM_MPC "Performs MPC simulation using RaiSim" OFF) 
if (SIM_MPC)
  find_package(raisim CONFIG REQUIRED)
  find_package(raisimOgre CONFIG REQUIRED)

  add_executable(
    mpc_standing
    mpc_standing.cpp
  )
  target_link_libraries(
    mpc_standing
    PRIVATE
    idocp::idocp
    raisim::raisim
    raisim::raisimOgre
  )
  target_include_directories(
    mpc_standing
    PRIVATE
    ${IDOCP_INCLUDE_DIR}
  )

  add_executable(
    mpc_trotting
    mpc_trotting.cpp
  )
  target_link_libraries(
    mpc_trotting
    PRIVATE
    idocp::idocp
    raisim::raisim
    raisim::raisimOgre
  )
  target_include_directories(
    mpc_trotting
    PRIVATE
    ${IDOCP_INCLUDE_DIR}
  )
endif()