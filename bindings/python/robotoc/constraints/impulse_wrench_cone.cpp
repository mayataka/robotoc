#include <pybind11/pybind11.h>

#include "robotoc/constraints/impulse_wrench_cone.hpp"
#include "robotoc/utils/pybind11_macros.hpp"


namespace robotoc {
namespace python {

namespace py = pybind11;

PYBIND11_MODULE(impulse_wrench_cone, m) {
  py::class_<ImpulseWrenchCone, ImpulseConstraintComponentBase, 
             std::shared_ptr<ImpulseWrenchCone>>(m, "ImpulseWrenchCone")
    .def(py::init<const Robot&, const double, const double>(),
          py::arg("robot"), py::arg("X"), py::arg("Y"))
    .def("set_rectangular", &ImpulseWrenchCone::setRectangular,
          py::arg("X"), py::arg("Y"))
    DEFINE_ROBOTOC_PYBIND11_CLASS_CLONE(ImpulseWrenchCone);
}

} // namespace python
} // namespace robotoc