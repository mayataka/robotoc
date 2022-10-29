cmake_minimum_required(VERSION 3.1)
project(robotoc_examples_iiwa14 CXX)

set(CMAKE_CXX_STANDARD 11)

find_package(robotoc REQUIRED)

macro(add_benchmark BENCHMARK)
  add_executable(
    ${BENCHMARK} 
    ${BENCHMARK}.cpp
  )
  target_link_libraries(
    ${BENCHMARK} 
    PRIVATE
    robotoc::robotoc
  )
endmacro()

macro(add_example EXAMPLE)
  add_executable(
    ${EXAMPLE} 
    ${EXAMPLE}.cpp
  )
  target_link_libraries(
    ${EXAMPLE} 
    PRIVATE
    robotoc::robotoc
  )
endmacro()


add_benchmark(unconstr_ocp_benchmark)
add_benchmark(unconstr_parnmpc_benchmark)
add_benchmark(ocp_benchmark)

add_example(config_space_ocp)
add_example(task_space_ocp)