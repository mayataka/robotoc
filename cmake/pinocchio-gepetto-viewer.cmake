include(FetchContent)
FetchContent_Declare(
  pinocchio-gepetto-viewer
  GIT_REPOSITORY https://github.com/mayataka/pinocchio-gepetto-viewer.git
)
FetchContent_Populate(pinocchio-gepetto-viewer)
FetchContent_MakeAvailable(pinocchio-gepetto-viewer)
set(external_pinocchio-gepetto-viewer_DIR ${PROJECT_BINARY_DIR}/_deps/pinocchio-gepetto-viewer-src)
configure_file(
  ${external_pinocchio-gepetto-viewer_DIR}/include/pinocchio/gepetto/viewer.hpp
  ${PROJECT_SOURCE_DIR}/include/idocp/third-party/pinocchio/gepetto/viewer.hpp
  COPYONLY
)
configure_file(
  ${external_pinocchio-gepetto-viewer_DIR}/include/pinocchio/gepetto/viewer.hxx
  ${PROJECT_SOURCE_DIR}/include/idocp/third-party/pinocchio/gepetto/viewer.hxx
  COPYONLY
)
configure_file(
  ${external_pinocchio-gepetto-viewer_DIR}/src/viewer.cc
  ${PROJECT_SOURCE_DIR}/src/third-party/viewer.cpp
  COPYONLY
)