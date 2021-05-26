pybind11_add_idocp_module(ocp_solver)
pybind11_add_idocp_module(unconstr_ocp_solver)
pybind11_add_idocp_module(unconstr_parnmpc_solver)

configure_file(
  __init__.py
  ${CMAKE_BINARY_DIR}/bindings/python/idocp/solver/__init__.py
  COPYONLY
)

# install(
#   FILES __init__.py
#   RUNTIME DESTINATION ${CMAKE_BINARY_DIR}/bindings/python/idocp/cost
# )
#   DESTINATION ${CMAKE_INSTALL_LIBDIR}/${PROJECT_NAME}
#   RUNTIME DESTINATION ${CMAKE_INSTALL_LIBDIR}/${PROJECT_NAME}
# )