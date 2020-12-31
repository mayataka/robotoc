cmake_minimum_required(VERSION 3.1)
project(idocp_examples_anymal CXX)

set(CMAKE_CXX_STANDARD 11)

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
    ocp_trotting
    ocp_trotting.cpp
)
target_link_libraries(
    ocp_trotting
    PRIVATE
    idocp::idocp
)
target_include_directories(
    ocp_trotting
    PRIVATE
    ${IDOCP_INCLUDE_DIR}
)