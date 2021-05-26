pybind11_add_idocp_module(cost_function_component_base)
pybind11_add_idocp_module(configuration_space_cost)
pybind11_add_idocp_module(task_space_3d_cost)
pybind11_add_idocp_module(task_space_6d_cost)
pybind11_add_idocp_module(com_cost)
pybind11_add_idocp_module(contact_force_cost)
pybind11_add_idocp_module(time_varying_configuration_ref_base)
pybind11_add_idocp_module(time_varying_configuration_space_cost)

configure_file(
  __init__.py
  ${CMAKE_BINARY_DIR}/bindings/python/idocp/cost/__init__.py
  COPYONLY
)

# install(
#   FILES __init__.py
#   RUNTIME DESTINATION ${CMAKE_BINARY_DIR}/bindings/python/idocp/cost
# )
#   DESTINATION ${CMAKE_INSTALL_LIBDIR}/${PROJECT_NAME}
#   RUNTIME DESTINATION ${CMAKE_INSTALL_LIBDIR}/${PROJECT_NAME}
# )