# add_executable(
#     joint_position_lower_limits_test
#     joint_position_lower_limits_test.cpp
# )
# target_include_directories(
#     joint_position_lower_limits_test
#     PRIVATE
#     ${GTEST_INCLUDE_PATH}
#     ${GMOCK_INCLUDE_PATH}
#     ${GTEST_LIBRARY_PATH}
#     ${EIGEN3_INCLUDE_DIR}
#     ${IDOCP_INCLUDE_DIR}
#     ${Boost_INCLUDE_DIRS}
#     ${PINOCCHIO_INCLUDEDIR}
# )
# target_link_libraries(
#     joint_position_lower_limits_test
#     PRIVATE
#     pthread
#     GTest::GTest
#     GTest::GMock
#     idocp_constraints
#     idocp_robot
# )
# add_dependencies(
#     joint_position_lower_limits_test
#     googletest
# )
# add_test(
#     NAME joint_position_lower_limits_test
#     COMMAND $<TARGET_FILE:joint_position_lower_limits_test>
# )

add_executable(
    joint_position_upper_limits_test
    joint_position_upper_limits_test.cpp
)
target_include_directories(
    joint_position_upper_limits_test
    PRIVATE
    ${GTEST_INCLUDE_PATH}
    ${GMOCK_INCLUDE_PATH}
    ${GTEST_LIBRARY_PATH}
    ${EIGEN3_INCLUDE_DIR}
    ${IDOCP_INCLUDE_DIR}
    ${Boost_INCLUDE_DIRS}
    ${PINOCCHIO_INCLUDEDIR}
)
target_link_libraries(
    joint_position_upper_limits_test
    PRIVATE
    pthread
    GTest::GTest
    GTest::GMock
    idocp_constraints
    idocp_robot
)
add_dependencies(
    joint_position_upper_limits_test
    googletest
)
add_test(
    NAME joint_position_upper_limits_test
    COMMAND $<TARGET_FILE:joint_position_upper_limits_test>
)

# add_executable(
#     joint_velocity_lower_limits_test
#     joint_velocity_lower_limits_test.cpp
# )
# target_include_directories(
#     joint_velocity_lower_limits_test
#     PRIVATE
#     ${GTEST_INCLUDE_PATH}
#     ${GMOCK_INCLUDE_PATH}
#     ${GTEST_LIBRARY_PATH}
#     ${EIGEN3_INCLUDE_DIR}
#     ${IDOCP_INCLUDE_DIR}
#     ${Boost_INCLUDE_DIRS}
#     ${PINOCCHIO_INCLUDEDIR}
# )
# target_link_libraries(
#     joint_velocity_lower_limits_test
#     PRIVATE
#     pthread
#     GTest::GTest
#     GTest::GMock
#     idocp_constraints
#     idocp_robot
# )
# add_dependencies(
#     joint_velocity_lower_limits_test
#     googletest
# )
# add_test(
#     NAME joint_velocity_lower_limits_test
#     COMMAND $<TARGET_FILE:joint_velocity_lower_limits_test>
# )

# add_executable(
#     joint_velocity_upper_limits_test
#     joint_velocity_upper_limits_test.cpp
# )
# target_include_directories(
#     joint_velocity_upper_limits_test
#     PRIVATE
#     ${GTEST_INCLUDE_PATH}
#     ${GMOCK_INCLUDE_PATH}
#     ${GTEST_LIBRARY_PATH}
#     ${EIGEN3_INCLUDE_DIR}
#     ${IDOCP_INCLUDE_DIR}
#     ${Boost_INCLUDE_DIRS}
#     ${PINOCCHIO_INCLUDEDIR}
# )
# target_link_libraries(
#     joint_velocity_upper_limits_test
#     PRIVATE
#     pthread
#     GTest::GTest
#     GTest::GMock
#     idocp_constraints
#     idocp_robot
# )
# add_dependencies(
#     joint_velocity_upper_limits_test
#     googletest
# )
# add_test(
#     NAME joint_velocity_upper_limits_test
#     COMMAND $<TARGET_FILE:joint_velocity_upper_limits_test>
# )