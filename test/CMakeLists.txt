# enable google test
include(cmake/gtest.cmake)

# build test helper 
add_subdirectory(test_helper)

# macto for testing
macro(add_robotoc_test TESTNAME)
  add_executable(
    ${TESTNAME} 
    ${TESTNAME}.cpp
  )
  target_include_directories(
    ${TESTNAME} 
    PRIVATE
    ${GTEST_INCLUDE_PATH}
    ${GMOCK_INCLUDE_PATH}
    ${GTEST_LIBRARY_PATH}
    ${PROJECT_SOURCE_DIR}/include
    ${PROJECT_SOURCE_DIR}/test/test_helper
    ${PINOCCHIO_INCLUDE_DIRS}
  )
  target_link_libraries(
    ${TESTNAME} 
    PRIVATE
    GTest::GTest
    GTest::GMock
    ${PROJECT_NAME}
    test_helper
  )
  add_dependencies(
    ${TESTNAME} 
    googletest
  )
  add_test(
    NAME ${TESTNAME} 
    COMMAND $<TARGET_FILE:${TESTNAME}>
  )
endmacro()

# copy urdf files
configure_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/urdf/iiwa14/iiwa14.urdf
  ${CMAKE_CURRENT_BINARY_DIR}/urdf/iiwa14/iiwa14.urdf
  COPYONLY
)
configure_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/urdf/anymal/anymal.urdf
  ${CMAKE_CURRENT_BINARY_DIR}/urdf/anymal/anymal.urdf
  COPYONLY
)
configure_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/urdf/icub/icub.urdf
  ${CMAKE_CURRENT_BINARY_DIR}/urdf/icub/icub.urdf
  COPYONLY
)

# add tests
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/robot)
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/core)
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/cost)
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/constraints)
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/dynamics)
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/ocp)
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/sto)
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/planner)
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/riccati)
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/unconstr)
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/parnmpc)
# add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/line_search)
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/solver)
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/mpc)