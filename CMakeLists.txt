cmake_minimum_required(VERSION 3.1)
project(idocp CXX)
set(CMAKE_CXX_STANDARD 11)

#############
## Options ##
#############
option(OPTIMIZE_FOR_NATIVE "Enable -march=native" OFF)
option(BUILD_VIEWER "Build trajectory viewer" OFF)
option(BUILD_TESTS "Build unit tests" OFF)

#################
## Build idocp ##
#################
if(NOT CMAKE_BUILD_TYPE MATCHES Release)
  message(STATUS "WARNING: CMAKE_BUILD_TYPE is NOT set to Release, which can decrease performance significantly.")
endif()
set(CMAKE_INSTALL_RPATH ${CMAKE_INSTALL_PREFIX}/lib)
# find Eigen
find_package(Eigen3 REQUIRED)
# find pinocchio
find_package(pinocchio REQUIRED)
# find OpenMP
find_package(OpenMP REQUIRED)
# build idocp 
add_library(
  ${PROJECT_NAME} 
  SHARED
  src/robot/robot.cpp
  src/robot/floating_base.cpp
  src/robot/point_contact.cpp
  src/robot/contact_status.cpp
  src/robot/impulse_status.cpp
  src/cost/cost_function.cpp
  src/cost/cost_function_data.cpp
  src/cost/cost_function_component_base.cpp
  src/cost/configuration_space_cost.cpp
  src/cost/task_space_3d_cost.cpp
  src/cost/task_space_6d_cost.cpp
  src/cost/com_cost.cpp
  src/cost/contact_force_cost.cpp
  src/cost/time_varying_configuration_space_cost.cpp
  src/cost/time_varying_task_space_3d_cost.cpp
  src/cost/time_varying_task_space_6d_cost.cpp
  src/cost/time_varying_com_cost.cpp
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
  src/constraints/impulse_friction_cone.cpp
  src/impulse/impulse_split_solution.cpp
  src/impulse/impulse_split_direction.cpp
  src/impulse/impulse_split_kkt_residual.cpp
  src/impulse/impulse_split_kkt_matrix.cpp
  src/impulse/impulse_state_equation.cpp
  src/impulse/impulse_dynamics_data.cpp
  src/impulse/impulse_dynamics.cpp
  src/impulse/impulse_split_ocp.cpp
  src/impulse/impulse_backward_riccati_recursion_factorizer.cpp
  src/impulse/impulse_split_riccati_factorizer.cpp
  src/hybrid/discrete_event.cpp
  src/hybrid/contact_sequence.cpp
  src/hybrid/ocp_discretizer.cpp
  src/hybrid/hybrid_container.cpp
  src/hybrid/collision_checker.cpp
  src/ocp/split_kkt_residual.cpp
  src/ocp/split_kkt_matrix.cpp
  src/ocp/split_solution.cpp
  src/ocp/split_direction.cpp
  src/ocp/state_equation.cpp
  src/ocp/contact_dynamics_data.cpp
  src/ocp/contact_dynamics.cpp
  src/ocp/switching_constraint.cpp
  src/ocp/split_state_constraint_jacobian.cpp
  src/ocp/switching_constraint.cpp
  src/ocp/split_ocp.cpp
  src/ocp/terminal_ocp.cpp
  src/ocp/state_constraint_jacobian.cpp
  src/ocp/ocp_linearizer.cpp
  src/ocp/split_riccati_factorization.cpp
  src/ocp/split_constrained_riccati_factorization.cpp
  src/ocp/backward_riccati_recursion_factorizer.cpp
  src/ocp/lqr_policy.cpp
  src/ocp/split_riccati_factorizer.cpp
  src/ocp/riccati_recursion_solver.cpp
  src/ocp/ocp_solver.cpp
  src/unocp/split_unkkt_residual.cpp
  src/unocp/split_unkkt_matrix.cpp
  src/unocp/split_unkkt_matrix_inverter.cpp
  src/unocp/unconstrained_container.cpp
  src/unocp/unconstrained_dynamics.cpp
  src/unocp/split_unocp.cpp
  src/unocp/backward_unriccati_recursion_factorizer.cpp
  src/unocp/split_unriccati_factorizer.cpp
  src/unocp/unriccati_recursion.cpp
  src/unocp/unocp_solver.cpp
  # src/unocp/split_unparnmpc.cpp
  # src/unocp/terminal_unparnmpc.cpp
  # src/unocp/split_unbackward_correction.cpp
  # src/unocp/unbackward_correction.cpp
  # src/unocp/unparnmpc_solver.cpp
  src/line_search/line_search_filter.cpp
  src/line_search/line_search.cpp
  src/line_search/unline_search.cpp
  src/utils/derivative_checker.cpp
  src/utils/joint_constraints_factory.cpp
  src/utils/ocp_benchmarker.cpp
  src/utils/logger.cpp
)
target_link_libraries(
  ${PROJECT_NAME} 
  PUBLIC
  ${PINOCCHIO_LIBRARIES}
  PRIVATE
  ${OpenMP_CXX_FLAGS}
)
target_include_directories(
  ${PROJECT_NAME} 
  PUBLIC
  ${EIGEN3_INCLUDE_DIR}
  ${PINOCCHIO_INCLUDE_DIRS}
  $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:include>
)
target_compile_options(
  ${PROJECT_NAME} 
  PRIVATE
  ${OpenMP_CXX_FLAGS}
)

##################
## Build viewer ##
##################
if (BUILD_VIEWER) 
  find_package(PkgConfig)
  pkg_check_modules(pinocchio-gepetto-viewer REQUIRED pinocchio-gepetto-viewer)
  target_sources(
    ${PROJECT_NAME} 
    PRIVATE
    src/utils/trajectory_viewer.cpp
  )
  target_link_libraries(
    ${PROJECT_NAME} 
    PUBLIC
    ${pinocchio-gepetto-viewer_LINK_LIBRARIES}
  )
  target_include_directories(
    ${PROJECT_NAME} 
    PUBLIC
    ${pinocchio-gepetto-viewer_INCLUDE_DIRS}
  )
endif() # end if (BUILD_VIEWER)

#############
## Testing ##
#############
# Add testing directory
if (BUILD_TESTS)
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
endif() # end if (BUILD_TESTS)

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