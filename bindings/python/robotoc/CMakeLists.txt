set(ROBOTOC_PYTHON_VERSION ${PYTHON_VERSION_MAJOR}.${PYTHON_VERSION_MINOR})
set(ROBOTOC_PYTHON_BINDINGS_LIBDIR ${CMAKE_INSTALL_LIBDIR}/python${ROBOTOC_PYTHON_VERSION}/site-packages/robotoc)

macro(pybind11_add_robotoc_module MODULE)
  pybind11_add_module(
    ${MODULE} 
    SHARED 
    ${MODULE}.cpp
  )
  target_link_libraries(
    ${MODULE} 
    PRIVATE
    robotoc::robotoc
  )
endmacro()


macro(install_robotoc_pybind_module CURRENT_MODULE_DIR)
  file(GLOB PYTHON_BINDINGS_${CURRENT_MODULE_DIR} ${CMAKE_CURRENT_BINARY_DIR}/*.cpython*)
  file(GLOB PYTHON_FILES_${CURRENT_MODULE_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/*.py)
  install(
    FILES ${PYTHON_FILES_${CURRENT_MODULE_DIR}} ${PYTHON_BINDINGS_${CURRENT_MODULE_DIR}} 
    DESTINATION ${ROBOTOC_PYTHON_BINDINGS_LIBDIR}/${CURRENT_MODULE_DIR}
  )
endmacro()

add_subdirectory(robot)
add_subdirectory(core)
add_subdirectory(cost)
add_subdirectory(constraints)
# add_subdirectory(dynamics)
add_subdirectory(ocp)
add_subdirectory(hybrid)
add_subdirectory(unconstr)
add_subdirectory(riccati)
add_subdirectory(line_search)
add_subdirectory(solver)
add_subdirectory(mpc)
add_subdirectory(utils)

install(
  FILES __init__.py 
  DESTINATION ${ROBOTOC_PYTHON_BINDINGS_LIBDIR}/
)