cmake_minimum_required(VERSION 3.5)
project(cranex7_examples CXX)

set(CMAKE_CXX_STANDARD 11)

find_package(idocp REQUIRED)

add_executable(
    cranex7_task_space_3d_ocp
    cranex7_task_space_3d_ocp.cpp
)
target_link_libraries(
    cranex7_task_space_3d_ocp
    PRIVATE
    idocp::idocp
)
target_include_directories(
    cranex7_task_space_3d_ocp
    PRIVATE
    ${IDOCP_INCLUDE_DIR}
)

add_executable(
    cranex7_task_space_6d_ocp
    cranex7_task_space_6d_ocp.cpp
)
target_link_libraries(
    cranex7_task_space_6d_ocp
    PRIVATE
    idocp::idocp
)
target_include_directories(
    cranex7_task_space_6d_ocp
    PRIVATE
    ${IDOCP_INCLUDE_DIR}
)

add_executable(
    cranex7_task_space_3d_parnmpc
    cranex7_task_space_3d_parnmpc.cpp
)
target_link_libraries(
    cranex7_task_space_3d_parnmpc
    PRIVATE
    idocp::idocp
)
target_include_directories(
    cranex7_task_space_3d_parnmpc
    PRIVATE
    ${IDOCP_INCLUDE_DIR}
)

add_executable(
    cranex7_task_space_6d_parnmpc
    cranex7_task_space_6d_parnmpc.cpp
)
target_link_libraries(
    cranex7_task_space_6d_parnmpc
    PRIVATE
    idocp::idocp
)
target_include_directories(
    cranex7_task_space_6d_parnmpc
    PRIVATE
    ${IDOCP_INCLUDE_DIR}
)