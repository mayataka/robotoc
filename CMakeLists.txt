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
# find Boost
find_package(Boost COMPONENTS system filesystem REQUIRED)
# find pinocchio
find_package(pinocchio REQUIRED)
# find OpenMP
find_package(OpenMP REQUIRED)
if(OpenMP_FOUND)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()
# build idocp 
set(IDOCP_SOURCE_FILES ${IDOCP_SOURCE_FILES} 
  src/robot/robot.cpp
  src/robot/floating_base.cpp
  src/robot/point_contact.cpp
  src/robot/contact_status.cpp
  src/robot/impulse_status.cpp
  src/cost/cost_function.cpp
  src/cost/cost_function_data.cpp
  src/cost/cost_function_component_base.cpp
  src/cost/configuration_space_cost.cpp
  src/cost/time_varying_configuration_space_cost.cpp
  src/cost/contact_force_cost.cpp
  src/cost/task_space_3d_cost.cpp
  src/cost/task_space_6d_cost.cpp
  src/cost/time_varying_task_space_3d_cost.cpp
  src/cost/trotting_configuration_cost.cpp
  src/cost/trotting_foot_step_cost.cpp
  src/constraints/pdipm.cpp
  src/constraints/constraints.cpp
  src/constraints/constraints_impl.cpp
  src/constraints/constraints_data.cpp
  src/constraints/constraint_component_base.cpp
  src/constraints/impulse_constraint_component_base.cpp
  src/constraints/constraint_component_data.cpp
  src/constraints/joint_position_upper_limit.cpp
  src/constraints/joint_position_lower_limit.cpp
  src/constraints/joint_velocity_upper_limit.cpp
  src/constraints/joint_velocity_lower_limit.cpp
  src/constraints/joint_acceleration_upper_limit.cpp
  src/constraints/joint_acceleration_lower_limit.cpp
  src/constraints/joint_torques_upper_limit.cpp
  src/constraints/joint_torques_lower_limit.cpp
  src/constraints/friction_cone.cpp
  src/constraints/contact_normal_force.cpp
  src/constraints/contact_distance.cpp
  src/constraints/impulse_friction_cone.cpp
  src/constraints/impulse_normal_force.cpp
  src/impulse/dynamic_schur_complement.cpp
  src/impulse/impulse_split_kkt_residual.cpp
  src/impulse/impulse_split_kkt_matrix.cpp
  src/impulse/impulse_split_solution.cpp
  src/impulse/impulse_split_direction.cpp
  src/impulse/impulse_state_equation.cpp
  src/impulse/impulse_dynamics_forward_euler_data.cpp
  src/impulse/impulse_dynamics_forward_euler.cpp
  src/impulse/impulse_split_ocp.cpp
  src/impulse/impulse_backward_riccati_recursion_factorizer.cpp
  src/impulse/impulse_forward_riccati_recursion_factorizer.cpp
  src/impulse/impulse_split_riccati_factorizer.cpp
  src/impulse/impulse_dynamics_backward_euler_data.cpp
  src/impulse/impulse_dynamics_backward_euler.cpp
  src/hybrid/discrete_event.cpp
  src/hybrid/contact_sequence.cpp
  src/hybrid/ocp_discretizer.cpp
  src/hybrid/hybrid_container.cpp
  src/hybrid/collision_checker.cpp
  src/ocp/schur_complement.cpp
  src/ocp/split_kkt_residual.cpp
  src/ocp/split_kkt_matrix.cpp
  src/ocp/split_solution.cpp
  src/ocp/split_direction.cpp
  src/ocp/state_equation.cpp
  src/ocp/contact_dynamics_data.cpp
  src/ocp/contact_dynamics.cpp
  src/ocp/split_ocp.cpp
  src/ocp/terminal_ocp.cpp
  src/ocp/ocp_linearizer.cpp
  src/ocp/split_riccati_factorization.cpp
  src/ocp/backward_riccati_recursion_factorizer.cpp
  src/ocp/forward_riccati_recursion_factorizer.cpp
  src/ocp/lqr_state_feedback_policy.cpp
  src/ocp/split_riccati_factorizer.cpp
  src/ocp/riccati_recursion.cpp
  src/ocp/split_state_constraint_riccati_factorization.cpp
  src/ocp/state_constraint_riccati_factorization.cpp
  src/ocp/state_constraint_riccati_lp_factorizer.cpp
  src/ocp/state_constraint_riccati_factorizer.cpp
  src/ocp/riccati_direction_calculator.cpp
  src/ocp/riccati_solver.cpp
  src/ocp/ocp_solver.cpp
  src/utils/joint_constraints_factory.cpp
  src/utils/ocp_benchmarker.cpp
  src/utils/simulation_data_saver.cpp
)
if (BUILD_WITH_RAISIM) 
  # find raisim and raisimOgre
  find_package(raisim CONFIG REQUIRED)
  find_package(raisimOgre CONFIG REQUIRED)
  list(APPEND IDOCP_SOURCE_FILES ${IDOCP_SOURCE_FILES}
    src/utils/helper.cpp
    src/utils/csv_to_eigen.cpp
    src/utils/trajectory_viewer.cpp
    # src/utils/quadruped_simulator.cpp
  )
endif() # end if (BUILD_WITH_RAISIM)
add_library(
  ${PROJECT_NAME} 
  SHARED
  ${IDOCP_SOURCE_FILES}
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
if (BUILD_WITH_RAISIM) 
  target_link_libraries(
    ${PROJECT_NAME} 
    PUBLIC 
    ${raisim_LIBRARIES}
    ${raisimOgre_LIBRARIES}
  )
  target_include_directories(
    ${PROJECT_NAME} 
    PUBLIC
    ${raisim_INCLUDE_DIRS}
    ${raisimOgre_INCLUDE_DIRS}
  )
endif() # end if (BUILD_WITH_RAISIM)

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