#include <pybind11/pybind11.h>

#include "robotoc/constraints/impact_wrench_cone.hpp"
#include "robotoc/utils/pybind11_macros.hpp"


namespace robotoc {
namespace python {

namespace py = pybind11;

PYBIND11_MODULE(impact_wrench_cone, m) {
  py::class_<ImpactWrenchCone, ImpactConstraintComponentBase, 
             std::shared_ptr<ImpactWrenchCone>>(m, "ImpactWrenchCone")
    .def(py::init<const Robot&, const double, const double>(),
          py::arg("robot"), py::arg("X"), py::arg("Y"))
    .def("set_rectangular", &ImpactWrenchCone::setRectangular,
          py::arg("X"), py::arg("Y"))
    DEFINE_ROBOTOC_PYBIND11_CLASS_CLONE(ImpactWrenchCone);
}

} // namespace python
} // namespace robotoc