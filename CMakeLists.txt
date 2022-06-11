cmake_minimum_required(VERSION 3.11)
project(robotoc CXX)

if (NOT "${CMAKE_CXX_STANDARD}")
  set(CMAKE_CXX_STANDARD 11)
else()
  if ("${CMAKE_CXX_STANDARD}" LESS 11)
    set(CMAKE_CXX_STANDARD 11)
  endif()
endif()

#############
## Options ##
#############
option(OPTIMIZE_FOR_NATIVE "Enable -march=native" OFF)
option(BUILD_VIEWER "Build trajectory viewer" OFF)
option(BUILD_TESTS "Build unit tests" OFF)
option(BUILD_PYTHON_INTERFACE "Build Python interface" ON)

###################
## Build robotoc ##
###################
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
# build robotoc 
file(GLOB_RECURSE ${PROJECT_NAME}_SOURCES src/*.cpp)
file(GLOB_RECURSE ${PROJECT_NAME}_HEADERS include/${PROJECT_NAME}/*.h*)
list(FILTER ${PROJECT_NAME}_SOURCES EXCLUDE REGEX ".trajectory_viewer.")
list(FILTER ${PROJECT_NAME}_HEADERS EXCLUDE REGEX ".trajectory_viewer.")
add_library(
  ${PROJECT_NAME} 
  SHARED
  ${${PROJECT_NAME}_SOURCES} 
  ${${PROJECT_NAME}_HEADERS}
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
if (OPTIMIZE_FOR_NATIVE)
  target_compile_options(
    ${PROJECT_NAME} 
    PUBLIC
    -march=native
  )
endif()

#############
## Mac OSX ##
#############
# In Mac OSX, we have to explicitly include Boost to avoid compilation errors.
if (${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
    set(MACOSX TRUE)
endif()
if (MACOSX)
  find_package(Boost REQUIRED)
  target_include_directories(
    ${PROJECT_NAME} 
    PUBLIC
    ${Boost_INCLUDE_DIRS}
  )
endif()

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

##############
## Bindings ##
##############
if (BUILD_PYTHON_INTERFACE)
  add_subdirectory(bindings/python)
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
set(CONFIG_PACKAGE_LOCATION "share/${PROJECT_NAME}/cmake")
install(
  EXPORT ${PROJECT_NAME}-config
  FILE ${PROJECT_NAME}-config.cmake
  NAMESPACE ${PROJECT_NAME}::
  DESTINATION ${CONFIG_PACKAGE_LOCATION}
)
# Install alias
add_library(
  ${PROJECT_NAME}::${PROJECT_NAME} 
  ALIAS ${PROJECT_NAME}
)