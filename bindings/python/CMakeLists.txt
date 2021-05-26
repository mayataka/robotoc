include(FetchContent)
FetchContent_GetProperties(pybind11)
if(NOT pybind11_POPULATED)
  FetchContent_Populate(
    pybind11
    GIT_REPOSITORY https://github.com/pybind/pybind11.git 
    GIT_TAG        master
    SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/pybind11
    SUBBUILD_DIR ${CMAKE_CURRENT_SOURCE_DIR}/pybind11-build)
endif()

add_subdirectory(pybind11)
add_subdirectory(idocp)