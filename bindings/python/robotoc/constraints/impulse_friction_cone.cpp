#include <pybind11/pybind11.h>

#include "robotoc/constraints/impulse_friction_cone.hpp"
#include "robotoc/utils/pybind11_macros.hpp"


namespace robotoc {
namespace python {

namespace py = pybind11;

PYBIND11_MODULE(impulse_friction_cone, m) {
  py::class_<ImpulseFrictionCone, ImpulseConstraintComponentBase, 
             std::shared_ptr<ImpulseFrictionCone>>(m, "ImpulseFrictionCone")
    .def(py::init<const Robot&>(),
          py::arg("robot"))
    DEFINE_ROBOTOC_PYBIND11_CLASS_CLONE(ImpulseFrictionCone);
}

} // namespace python
} // namespace robotoc