include(ExternalProject)

# find pthread for googletest 
find_package(Threads REQUIRED)

SET_DIRECTORY_PROPERTIES(PROPERTIES EP_PREFIX ${CMAKE_BINARY_DIR}/external)
externalproject_add(
    googletest
    URL https://github.com/google/googletest/archive/release-1.10.0.zip
    UPDATE_COMMAND ""
    INSTALL_COMMAND ""
    )

externalproject_get_property(googletest source_dir)
set(GTEST_INCLUDE_PATH ${source_dir}/googletest/include)
set(GMOCK_INCLUDE_PATH ${source_dir}/googlemock/include)

externalproject_get_property(googletest binary_dir)
set(GTEST_LIBRARY_PATH ${binary_dir}/lib/${CMAKE_FIND_LIBRARY_PREFIXES}gtest.a) # in Unix, libgtest.a
# set(GTEST_LIBRARY_PATH ${binary_dir}/googlemock/gtest/${CMAKE_FIND_LIBRARY_PREFIXES}gtest.a) # in Unix, libgtest.a
set(GTEST_LIBRARY GTest::GTest)
add_library(${GTEST_LIBRARY} UNKNOWN IMPORTED)
set_target_properties(${GTEST_LIBRARY} PROPERTIES
    IMPORTED_LOCATION ${GTEST_LIBRARY_PATH}
    INTERFACE_LINK_LIBRARIES Threads::Threads
    )
add_dependencies(${GTEST_LIBRARY} googletest)

set(GMOCK_LIBRARY_PATH ${binary_dir}/lib/${CMAKE_FIND_LIBRARY_PREFIXES}gmock.a) # in Unix, libgmock.a
# set(GMOCK_LIBRARY_PATH ${binary_dir}/googlemock/${CMAKE_FIND_LIBRARY_PREFIXES}gmock.a) # in Unix, libgmock.a
set(GMOCK_LIBRARY GTest::GMock)
add_library(${GMOCK_LIBRARY} UNKNOWN IMPORTED)
set_target_properties(${GMOCK_LIBRARY} PROPERTIES
    IMPORTED_LOCATION ${GMOCK_LIBRARY_PATH}
    INTERFACE_LINK_LIBRARIES Threads::Threads
    )
add_dependencies(${GMOCK_LIBRARY} googletest)