pybind11_add_idocp_module(split_solution)

configure_file(
  __init__.py
  ${CMAKE_BINARY_DIR}/bindings/python/idocp/ocp/__init__.py
  COPYONLY
)


# install(
#   FILES __init__.py
#   RUNTIME DESTINATION ${CMAKE_BINARY_DIR}/bindings/python/idocp/cost
# )
#   DESTINATION ${CMAKE_INSTALL_LIBDIR}/${PROJECT_NAME}
#   RUNTIME DESTINATION ${CMAKE_INSTALL_LIBDIR}/${PROJECT_NAME}
# )