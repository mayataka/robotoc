cmake_minimum_required(VERSION 3.1)
project(idocp CXX)
set(CMAKE_CXX_STANDARD 14)

###########
## Build ##
###########
if(NOT CMAKE_BUILD_TYPE MATCHES Release)
  message(STATUS "WARNING: CMAKE_BUILD_TYPE not set to Release, which can decrease performance significantly")
endif()
set(CMAKE_INSTALL_RPATH ${CMAKE_INSTALL_PREFIX}/lib)
# find Eigen
find_package(Eigen3 REQUIRED)
include_directories()
# find Boost
find_package(Boost COMPONENTS system filesystem REQUIRED)
# find pinocchio
find_package(PkgConfig)
pkg_check_modules(PINOCCHIO REQUIRED pinocchio)
add_definitions(
  -DPINOCCHIO_URDFDOM_TYPEDEF_SHARED_PTR
  -DPINOCCHIO_URDFDOM_USE_STD_SHARED_PTR
  -DPINOCCHIO_WITH_URDFDOM
)
link_directories(
  ${PINOCCHIO_LIBDIR}
)
# find OpenMP
find_package(OpenMP REQUIRED)
if(OpenMP_FOUND)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()
# build idocp 
add_library(
  ${PROJECT_NAME} 
  SHARED
  src/robot/robot.cpp
  src/robot/floating_base.cpp
  src/robot/point_contact.cpp
  # src/cost/cost_function.cpp
  # src/cost/cost_function_data.cpp
  # src/cost/joint_space_cost.cpp
  # src/cost/contact_cost.cpp
  # src/cost/task_space_3d_cost.cpp
  # src/cost/task_space_6d_cost.cpp
  # src/constraints/pdipm_func.cpp
  # src/constraints/constraints.cpp
  # src/constraints/constraint_component_base.cpp
  # src/constraints/constraint_component_data.cpp
  # src/constraints/joint_position_upper_limit.cpp
  # src/constraints/joint_position_lower_limit.cpp
  # src/constraints/joint_velocity_upper_limit.cpp
  # src/constraints/joint_velocity_lower_limit.cpp
  # src/constraints/joint_acceleration_upper_limit.cpp
  # src/constraints/joint_acceleration_lower_limit.cpp
  # src/constraints/joint_torques_upper_limit.cpp
  # src/constraints/joint_torques_lower_limit.cpp
  # src/ocp/kkt_residual.cpp
  # src/ocp/kkt_matrix.cpp
  # src/ocp/split_solution.cpp
  # src/ocp/split_direction.cpp
  # src/ocp/state_equation.cpp
  # src/ocp/robot_dynamics.cpp
  # src/ocp/line_search_filter.cpp
  # src/ocp/riccati_factorization.cpp
  # src/ocp/riccati_gain.cpp
  # src/ocp/riccati_matrix_factorizer.cpp
  # src/ocp/riccati_matrix_inverter.cpp
  # src/ocp/terminal_ocp.cpp
  # src/ocp/split_ocp.cpp
  # src/ocp/ocp.cpp
  # src/ocp/split_parnmpc.cpp
  # src/ocp/parnmpc.cpp
  # src/ocp/mpc.cpp
  # src/utils/joint_constraints_factory.cpp
  # src/utils/ocp_benchmarker.cpp
  # src/utils/runge_kutta.cpp
  # src/utils/simulation_data_saver.cpp
  # src/utils/simulator.cpp
  # src/utils/raisim_wrapper.cpp
)
target_link_libraries(
  ${PROJECT_NAME} 
  PUBLIC 
  pthread
  ${Boost_LIBRARIES}
  ${PINOCCHIO_LIBRARIES}
  urdfdom_model
)
target_include_directories(
  ${PROJECT_NAME} 
  PUBLIC
  ${EIGEN3_INCLUDE_DIR}
  ${Boost_INCLUDE_DIR}
  ${PINOCCHIO_INCLUDE_DIRS}
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:include>
)
set_target_properties(
  ${PROJECT_NAME} 
  PROPERTIES
  INTERFACE_LINK_DIRECTORIES $<INSTALL_INTERFACE:${PINOCCHIO_LIBDIR}>
)

#############
## Testing ##
#############
# Add testing directory
if (TESTING)
  enable_testing()
  # coverage options
  set(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} --coverage")
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} --coverage")
  add_library(coverage_config INTERFACE)
  if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    target_compile_options(coverage_config INTERFACE
      -O0        
      -g         
      --coverage 
    )
    if(CMAKE_VERSION VERSION_GREATER_EQUAL 3.13)
      target_link_options(coverage_config INTERFACE --coverage)
    else()
      target_link_libraries(coverage_config INTERFACE --coverage)
    endif()
  endif()
  add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/test)
endif() # end if (TESTING)

#############
## Install ##
#############
include(GNUInstallDirs)
# Install lib files
install(
  TARGETS ${PROJECT_NAME}
  EXPORT ${PROJECT_NAME}-config
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  RUNTIME DESTINATION ${CMAKE_INSTALL_LIBDIR}/${PROJECT_NAME}
)
# Install header files
install(
  DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/include/${PROJECT_NAME}/
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/${PROJECT_NAME}
)
# Install config file 
set(config_package_location "share/${PROJECT_NAME}/cmake")
install(
  EXPORT ${PROJECT_NAME}-config
  FILE ${PROJECT_NAME}-config.cmake
  NAMESPACE ${PROJECT_NAME}::
  DESTINATION ${config_package_location}
)
# Install alias
add_library(
  ${PROJECT_NAME}::${PROJECT_NAME} 
  ALIAS ${PROJECT_NAME}
)