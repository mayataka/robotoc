cmake_minimum_required(VERSION 3.5)
project(iiwa14 CXX)

set(CMAKE_CXX_STANDARD 14)

find_package(idocp REQUIRED)

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

add_executable(
    mpc_simulation
    mpc_simulation.cpp
)
target_link_libraries(
    mpc_simulation
    PRIVATE
    idocp::idocp
)
target_include_directories(
    mpc_simulation
    PRIVATE
    ${IDOCP_INCLUDE_DIR}
)