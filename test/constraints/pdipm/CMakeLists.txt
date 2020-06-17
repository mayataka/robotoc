add_executable(
    pdipm_func_test
    pdipm_func_test.cpp
    ${IDOCP_SRC_DIR}/constraints/pdipm/pdipm_func.cpp
)
target_include_directories(
    pdipm_func_test
    PRIVATE
    ${GTEST_INCLUDE_PATH}
    ${GMOCK_INCLUDE_PATH}
    ${GTEST_LIBRARY_PATH}
    ${IDOCP_INCLUDE_DIR}
)
target_link_libraries(
    pdipm_func_test
    PRIVATE
    pthread
    GTest::GTest
    GTest::GMock
)
add_dependencies(
    pdipm_func_test
    googletest
)
add_test(
    NAME pdipm_func_test
    COMMAND $<TARGET_FILE:pdipm_func_test>
)
