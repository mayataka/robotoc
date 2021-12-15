cmake_minimum_required(VERSION 3.1)
project(robotoc_examples_ANYmal CXX)

set(CMAKE_CXX_STANDARD 11)

find_package(robotoc REQUIRED)

option(ENABLE_VIEWER "Enable trajectory viewer" OFF)

macro(add_benchmark BENCHMARK)
  add_executable(
    ${BENCHMARK} 
    ${BENCHMARK}.cpp
  )
  target_include_directories(
    ${BENCHMARK} 
    PRIVATE
    ${ROBOTOC_INCLUDE_DIR}
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
  target_include_directories(
    ${EXAMPLE} 
    PRIVATE
  ${ROBOTOC_INCLUDE_DIR}
  )
  target_link_libraries(
    ${EXAMPLE} 
    PRIVATE
    robotoc::robotoc
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

add_example(trotting)
add_example(walking)
add_example(pacing)
add_example(bounding)
add_example(jumping)
add_example(running)

add_example(jumping_sto)