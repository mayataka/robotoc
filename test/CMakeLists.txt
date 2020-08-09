# enable google test
include(cmake/gtest.cmake)

# macto for testing
macro(add_test_macro TESTNAME)
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
    ${PINOCCHIO_INCLUDE_DIRS}
  )
  target_link_libraries(
    ${TESTNAME} 
    PRIVATE
    GTest::GTest
    GTest::GMock
    ${PROJECT_NAME}
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

# add tests
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/robot)
# add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/cost)
# add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/constraints)
# add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/manipulator)
# add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/ocp)