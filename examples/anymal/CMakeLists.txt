cmake_minimum_required(VERSION 3.1)
project(idocp_examples_ANYmal CXX)

set(CMAKE_CXX_STANDARD 11)

find_package(idocp REQUIRED)

macro(add_benchmark BENCHMARK)
  add_executable(
    ${BENCHMARK} 
    ${BENCHMARK}.cpp
  )
  target_include_directories(
    ${BENCHMARK} 
    PRIVATE
    ${IDOCP_INCLUDE_DIR}
  )
  target_link_libraries(
    ${BENCHMARK} 
    PRIVATE
    idocp::idocp
  )
endmacro()

add_library(
  anymal_costs  
  SHARED
  periodic_foot_track_ref.cpp
  periodic_com_ref.cpp
)
target_link_libraries(
  anymal_costs  
  PRIVATE
  idocp::idocp
)
target_include_directories(
  anymal_costs  
  PRIVATE
  ${IDOCP_INCLUDE_DIR}
)

macro(add_example EXAMPLE)
  add_executable(
    ${EXAMPLE} 
    ${EXAMPLE}.cpp
  )
  target_include_directories(
    ${EXAMPLE} 
    PRIVATE
  ${IDOCP_INCLUDE_DIR}
  )
  target_link_libraries(
    ${EXAMPLE} 
    PRIVATE
    idocp::idocp
    anymal_costs
  )
  if (ENABLE_VIEWER)
    target_compile_definitions(
      ${EXAMPLE} 
      PRIVATE
      ENABLE_VIEWER    
    )
  endif()
endmacro()

add_benchmark(ocp_benchmark)
add_benchmark(parnmpc_benchmark)
add_example(anymal_trotting)
# add_example(anymal_trotting_parnmpc)
# add_example(anymal_jumping)
# add_example(anymal_running)