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
link_directories(${PINOCCHIO_LIBDIR})

# set idocp directories
set(IDOCP_SRC_DIR ${CMAKE_CURRENT_SOURCE_DIR}/src)
set(IDOCP_INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/include/${PROJECT_NAME})

# build library
add_library(
  ${PROJECT_NAME} 
  SHARED
  src/robot/robot.cpp
  src/robot/floating_base.cpp
  src/robot/point_contact.cpp
  src/cost/joint_space_cost.cpp
  src/cost/contact_cost.cpp
  src/constraints/pdipm_func.cpp
  src/constraints/joint_space_constraints/joint_space_constraints.cpp
  src/constraints/joint_space_constraints/joint_variables_upper_limits.cpp
  src/constraints/joint_space_constraints/joint_variables_lower_limits.cpp
  src/manipulator/cost_function.cpp
  src/manipulator/constraints.cpp
  src/quadruped/cost_function.cpp
  src/quadruped/constraints.cpp
  src/ocp/riccati_matrix_factorizer.cpp
  src/ocp/riccati_matrix_inverter.cpp
  src/ocp/ocp_linearizer.cpp
  src/ocp/split_terminal_ocp.cpp
  src/ocp/split_ocp.cpp
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

#############
## Testing ##
#############
# Add testing directory
enable_testing()
if (TESTING)
  add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/test)
endif()

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