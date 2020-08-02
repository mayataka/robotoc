add_executable(
    riccati_matrix_factorizer_test
    riccati_matrix_factorizer_test.cpp
)
target_include_directories(
    riccati_matrix_factorizer_test
    PRIVATE
    ${GTEST_INCLUDE_PATH}
    ${GMOCK_INCLUDE_PATH}
    ${GTEST_LIBRARY_PATH}
    ${IDOCP_INCLUDE_DIR}
)
target_link_libraries(
    riccati_matrix_factorizer_test
    PRIVATE
    pthread
    GTest::GTest
    GTest::GMock
    idocp_robot 
    idocp_ocp
)
add_dependencies(
    riccati_matrix_factorizer_test
    googletest
)
add_test(
    NAME riccati_matrix_factorizer_test
    COMMAND $<TARGET_FILE:riccati_matrix_factorizer_test>
)

add_executable(
    riccati_matrix_inverter_test
    riccati_matrix_inverter_test.cpp
)
target_include_directories(
    riccati_matrix_inverter_test
    PRIVATE
    ${GTEST_INCLUDE_PATH}
    ${GMOCK_INCLUDE_PATH}
    ${GTEST_LIBRARY_PATH}
    ${IDOCP_INCLUDE_DIR}
)
target_link_libraries(
    riccati_matrix_inverter_test
    PRIVATE
    pthread
    GTest::GTest
    GTest::GMock
    idocp_robot 
    idocp_ocp
)
add_dependencies(
    riccati_matrix_inverter_test
    googletest
)
add_test(
    NAME riccati_matrix_inverter_test 
    COMMAND $<TARGET_FILE:riccati_matrix_inverter_test>
)

add_executable(
    split_terminal_ocp_test
    split_terminal_ocp_test.cpp
)
target_include_directories(
    split_terminal_ocp_test
    PRIVATE
    ${GTEST_INCLUDE_PATH}
    ${GMOCK_INCLUDE_PATH}
    ${GTEST_LIBRARY_PATH}
    ${IDOCP_INCLUDE_DIR}
)
target_link_libraries(
    split_terminal_ocp_test
    PRIVATE
    pthread
    GTest::GTest
    GTest::GMock
    idocp_robot 
    idocp_ocp
)
add_dependencies(
    split_terminal_ocp_test
    googletest
)
add_test(
    NAME split_terminal_ocp_test
    COMMAND $<TARGET_FILE:split_terminal_ocp_test>
)


add_executable(
    fixed_base_split_ocp_test
    fixed_base_split_ocp_test.cpp
)
target_include_directories(
    fixed_base_split_ocp_test
    PRIVATE
    ${GTEST_INCLUDE_PATH}
    ${GMOCK_INCLUDE_PATH}
    ${GTEST_LIBRARY_PATH}
    ${IDOCP_INCLUDE_DIR}
)
target_link_libraries(
    fixed_base_split_ocp_test
    PRIVATE
    pthread
    GTest::GTest
    GTest::GMock
    idocp_robot 
    idocp_ocp
)
add_dependencies(
    fixed_base_split_ocp_test
    googletest
)
add_test(
    NAME fixed_base_split_ocp_test
    COMMAND $<TARGET_FILE:fixed_base_split_ocp_test>
)