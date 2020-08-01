add_executable(
    floating_base_point_contact_test
    floating_base_point_contact_test.cpp
)
target_include_directories(
    floating_base_point_contact_test
    PRIVATE
    ${GTEST_INCLUDE_PATH}
    ${GMOCK_INCLUDE_PATH}
    ${GTEST_LIBRARY_PATH}
    ${IDOCP_INCLUDE_DIR}
)
target_link_libraries(
    floating_base_point_contact_test
    PRIVATE
    pthread
    GTest::GTest
    GTest::GMock
    idocp_robot
)
add_dependencies(
    floating_base_point_contact_test
    googletest
)
add_test(
    NAME floating_base_point_contact_test
    COMMAND $<TARGET_FILE:floating_base_point_contact_test>
)

add_executable(
    floating_base_floating_base_test
    floating_base_floating_base_test.cpp
)
target_include_directories(
    floating_base_floating_base_test
    PRIVATE
    ${GTEST_INCLUDE_PATH}
    ${GMOCK_INCLUDE_PATH}
    ${GTEST_LIBRARY_PATH}
    ${IDOCP_INCLUDE_DIR}
)
target_link_libraries(
    floating_base_floating_base_test
    PRIVATE
    pthread
    GTest::GTest
    GTest::GMock
    idocp_robot
)
add_dependencies(
    floating_base_floating_base_test
    googletest
)
add_test(
    NAME floating_base_floating_base_test
    COMMAND $<TARGET_FILE:floating_base_floating_base_test>
)

add_executable(
    floating_base_robot_test 
    floating_base_robot_test.cpp
)
target_include_directories(
    floating_base_robot_test 
    PRIVATE
    ${GTEST_INCLUDE_PATH}
    ${GMOCK_INCLUDE_PATH}
    ${GTEST_LIBRARY_PATH}
    ${IDOCP_INCLUDE_DIR}
)
target_link_libraries(
    floating_base_robot_test 
    PRIVATE
    pthread
    GTest::GTest
    GTest::GMock
    robot
)
add_dependencies(
    floating_base_robot_test 
    googletest
)
add_test(
    NAME floating_base_robot_test
    COMMAND $<TARGET_FILE:floating_base_robot_test>
)