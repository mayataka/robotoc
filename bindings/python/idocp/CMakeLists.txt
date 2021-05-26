macro(pybind11_add_idocp_module MODULE)
  pybind11_add_module(
    ${MODULE} 
    SHARED 
    ${MODULE}.cpp
  )
  target_link_libraries(
    ${MODULE} 
    PRIVATE
    idocp::idocp
  )
endmacro()

add_subdirectory(robot)
add_subdirectory(cost)
add_subdirectory(constraints)
add_subdirectory(ocp)
add_subdirectory(solver)

configure_file(
  __init__.py
  ${CMAKE_BINARY_DIR}/bindings/python/idocp/__init__.py
  COPYONLY
)