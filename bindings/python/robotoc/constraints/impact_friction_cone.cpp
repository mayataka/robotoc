#include <pybind11/pybind11.h>

#include "robotoc/constraints/impact_friction_cone.hpp"
#include "robotoc/utils/pybind11_macros.hpp"


namespace robotoc {
namespace python {

namespace py = pybind11;

PYBIND11_MODULE(impact_friction_cone, m) {
  py::class_<ImpactFrictionCone, ImpactConstraintComponentBase, 
             std::shared_ptr<ImpactFrictionCone>>(m, "ImpactFrictionCone")
    .def(py::init<const Robot&>(),
          py::arg("robot"))
    DEFINE_ROBOTOC_PYBIND11_CLASS_CLONE(ImpactFrictionCone);
}

} // namespace python
} // namespace robotoc