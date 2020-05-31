add_executable(
    passive_joints_test
    passive_joints_test.cpp
    ${SRC_DIR}/robot/passive_joints.cpp
)
target_include_directories(
    passive_joints_test
    PRIVATE
    ${GTEST_INCLUDE_PATH}
    ${GMOCK_INCLUDE_PATH}
    ${GTEST_LIBRARY_PATH}
    ${INCLUDE_DIR}
    ${Boost_INCLUDE_DIR}
    ${PINOCCHIO_INCLUDE_DIRS}
)
target_link_libraries(
    passive_joints_test
    PRIVATE
    pthread
    GTest::GTest
    GTest::GMock
    ${Boost_LIBRARIES}
    ${PINOCCHIO_LIBRARIES}
    urdfdom_model
)
add_dependencies(
    passive_joints_test
    googletest
)
add_test(
    NAME passive_joints_test
    COMMAND $<TARGET_FILE:passive_joints_test>
)


add_executable(
    point_contact_test
    point_contact_test.cpp
    ${SRC_DIR}/robot/point_contact.cpp
)
target_include_directories(
    point_contact_test
    PRIVATE
    ${GTEST_INCLUDE_PATH}
    ${GMOCK_INCLUDE_PATH}
    ${GTEST_LIBRARY_PATH}
    ${INCLUDE_DIR}
    ${Boost_INCLUDE_DIR}
    ${PINOCCHIO_INCLUDE_DIRS}
)
target_link_libraries(
    point_contact_test
    PRIVATE
    pthread
    GTest::GTest
    GTest::GMock
    ${Boost_LIBRARIES}
    ${PINOCCHIO_LIBRARIES}
    urdfdom_model
)
add_dependencies(
    point_contact_test
    googletest
)
add_test(
    NAME point_contact_test
    COMMAND $<TARGET_FILE:point_contact_test>
)


add_executable(
    fixed_base_robot_test 
    fixed_base_robot_test.cpp
    ${SRC_DIR}/robot/robot.cpp
    ${SRC_DIR}/robot/point_contact.cpp
    ${SRC_DIR}/robot/passive_joints.cpp
)
target_include_directories(
    fixed_base_robot_test 
    PRIVATE
    ${GTEST_INCLUDE_PATH}
    ${GMOCK_INCLUDE_PATH}
    ${GTEST_LIBRARY_PATH}
    ${INCLUDE_DIR}
    ${Boost_INCLUDE_DIR}
    ${PINOCCHIO_INCLUDE_DIRS}
)
target_link_libraries(
    fixed_base_robot_test 
    PRIVATE
    pthread
    GTest::GTest
    GTest::GMock
    ${Boost_LIBRARIES}
    ${PINOCCHIO_LIBRARIES}
    urdfdom_model
)
add_dependencies(
    fixed_base_robot_test 
    googletest
)
add_test(
    NAME fixed_base_robot_test
    COMMAND $<TARGET_FILE:fixed_base_robot_test>
)