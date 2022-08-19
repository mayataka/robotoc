#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "robotoc/constraints/friction_cone.hpp"


namespace robotoc {
namespace python {

namespace py = pybind11;

PYBIND11_MODULE(friction_cone, m) {
  py::class_<FrictionCone, ConstraintComponentBase, 
             std::shared_ptr<FrictionCone>>(m, "FrictionCone")
    .def(py::init<const Robot&>(),
          py::arg("robot"));
}

} // namespace python
} // namespace robotoc