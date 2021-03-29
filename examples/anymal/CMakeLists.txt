cmake_minimum_required(VERSION 3.1)
project(idocp_examples_ANYmal CXX)

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

# add_executable(
#     anymal_trotting
#     anymal_trotting.cpp
# )
# target_link_libraries(
#     anymal_trotting
#     PRIVATE
#     idocp::idocp
# )
# target_include_directories(
#     anymal_trotting
#     PRIVATE
#     ${IDOCP_INCLUDE_DIR}
# )
# if (ENABLE_VIEWER)
#     target_compile_definitions(
#         anymal_trotting
#         PRIVATE
#         ENABLE_VIEWER    
#     )
# endif()

# add_executable(
#     anymal_jumping
#     anymal_jumping.cpp
# )
# target_link_libraries(
#     anymal_jumping
#     PRIVATE
#     idocp::idocp
# )
# target_include_directories(
#     anymal_jumping
#     PRIVATE
#     ${IDOCP_INCLUDE_DIR}
# )
# if (ENABLE_VIEWER)
#     target_compile_definitions(
#         anymal_jumping
#         PRIVATE
#         ENABLE_VIEWER    
#     )
# endif()

# add_executable(
#     anymal_running
#     anymal_running.cpp
# )
# target_link_libraries(
#     anymal_running
#     PRIVATE
#     idocp::idocp
# )
# target_include_directories(
#     anymal_running
#     PRIVATE
#     ${IDOCP_INCLUDE_DIR}
# )
# if (ENABLE_VIEWER)
#     target_compile_definitions(
#         anymal_running
#         PRIVATE
#         ENABLE_VIEWER    
#     )
# endif()

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
if (ENABLE_VIEWER)
  target_compile_definitions(
    parnmpc_benchmark
    PRIVATE
    ENABLE_VIEWER    
  )
endif()

# add_executable(
#     anymal_trotting_parnmpc
#     anymal_trotting_parnmpc.cpp
# )
# target_link_libraries(
#     anymal_trotting_parnmpc
#     PRIVATE
#     idocp::idocp
# )
# target_include_directories(
#     anymal_trotting_parnmpc
#     PRIVATE
#     ${IDOCP_INCLUDE_DIR}
# )
# if (ENABLE_VIEWER)
#     target_compile_definitions(
#         anymal_trotting_parnmpc
#         PRIVATE
#         ENABLE_VIEWER    
#     )
# endif()