cmake_minimum_required(VERSION 3.1)
project(idocp_examples_iiwa14 CXX)

set(CMAKE_CXX_STANDARD 11)

find_package(idocp REQUIRED)

add_executable(
  unocp_benchmark
  unocp_benchmark.cpp
)
target_link_libraries(
  unocp_benchmark
  PRIVATE
  idocp::idocp
)
target_include_directories(
  unocp_benchmark
  PRIVATE
  ${IDOCP_INCLUDE_DIR}
)

add_executable(
  unparnmpc_benchmark
   unparnmpc_benchmark.cpp
)
target_link_libraries(
  unparnmpc_benchmark
  PRIVATE
  idocp::idocp
)
target_include_directories(
  unparnmpc_benchmark
  PRIVATE
  ${IDOCP_INCLUDE_DIR}
)

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
  parnmpc_benchmark
  parnmpc_benchmark.cpp
)
target_link_libraries(
  parnmpc_benchmark
  PRIVATE
  idocp::idocp
)
target_include_directories(
  parnmpc_benchmark
  PRIVATE
  ${IDOCP_INCLUDE_DIR}
)

add_executable(
  config_space_ocp
  config_space_ocp.cpp
)
target_link_libraries(
  config_space_ocp
  PRIVATE
  idocp::idocp
)
target_include_directories(
  config_space_ocp
  PRIVATE
  ${IDOCP_INCLUDE_DIR}
)
if (ENABLE_VIEWER)
  target_compile_definitions(
    config_space_ocp
    PRIVATE
    ENABLE_VIEWER    
  )
endif()

add_executable(
  task_space_ocp
  task_space_ocp.cpp
)
target_link_libraries(
  task_space_ocp
  PRIVATE
  idocp::idocp
)
target_include_directories(
  task_space_ocp
  PRIVATE
  ${IDOCP_INCLUDE_DIR}
)
if (ENABLE_VIEWER)
  target_compile_definitions(
    task_space_ocp
    PRIVATE
    ENABLE_VIEWER    
  )
endif()