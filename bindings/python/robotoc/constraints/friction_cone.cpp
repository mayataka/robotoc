#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "robotoc/constraints/friction_cone.hpp"
#include "robotoc/utils/pybind11_macros.hpp"


namespace robotoc {
namespace python {

namespace py = pybind11;

PYBIND11_MODULE(friction_cone, m) {
  py::class_<FrictionCone, ConstraintComponentBase, 
             std::shared_ptr<FrictionCone>>(m, "FrictionCone")
    .def(py::init<const Robot&>(),
          py::arg("robot"))
    DEFINE_ROBOTOC_PYBIND11_CLASS_CLONE(FrictionCone);
}

} // namespace python
} // namespace robotoc