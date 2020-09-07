cmake_minimum_required(VERSION 3.5)
project(anymal CXX)

set(CMAKE_CXX_STANDARD 14)

find_package(idocp REQUIRED)

find_package(PkgConfig)
pkg_check_modules(PINOCCHIO REQUIRED pinocchio)
link_directories(${PINOCCHIO_LIBDIR})

add_executable(
    ocp_benchmark
    ocp_benchmark.cpp
)
target_link_libraries(
    ocp_benchmark
    PRIVATE
    idocp::idocp
)
target_include_directories(
    ocp_benchmark
    PRIVATE
    ${IDOCP_INCLUDE_DIR}
)

option(SIM_MPC "Performs MPC simulation using RaiSim" OFF) 
set(SIM_MPC OFF)
if (SIM_MPC)
  find_package(raisim CONFIG REQUIRED)
  find_package(raisimOgre CONFIG REQUIRED)
  add_executable(
    mpc_simulation
    mpc_simulation.cpp
  )
  target_link_libraries(
    mpc_simulation
    PRIVATE
    idocp::idocp
    raisim::raisim
    raisim::raisimOgre
  )
  target_include_directories(
    mpc_simulation
    PRIVATE
    ${IDOCP_INCLUDE_DIR}
  )
endif()