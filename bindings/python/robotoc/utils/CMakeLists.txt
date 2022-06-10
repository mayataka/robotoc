pybind11_add_robotoc_module(openmp)
target_link_libraries(openmp PRIVATE ${OpenMP_CXX_FLAGS})

install_robotoc_pybind_module(utils)