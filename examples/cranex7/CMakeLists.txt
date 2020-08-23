cmake_minimum_required(VERSION 3.5)
project(crane7 CXX)

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