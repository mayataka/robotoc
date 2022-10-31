#include <pybind11/pybind11.h>

#include "robotoc/constraints/constraints.hpp"
#include "robotoc/utils/pybind11_macros.hpp"


namespace robotoc {
namespace python {

namespace py = pybind11;

using ConstraintComponentBasePtr = std::shared_ptr<ConstraintComponentBase>;
using ImpactConstraintComponentBasePtr = std::shared_ptr<ImpactConstraintComponentBase>;

PYBIND11_MODULE(constraints, m) {
  py::class_<Constraints, std::shared_ptr<Constraints>>(m, "Constraints")
    .def(py::init<const double, const double>(),
          py::arg("barrier_param")=1.0e-03, py::arg("fraction_to_boundary_rule")=0.995)
    .def("exist", &Constraints::exist,
          py::arg("name"))
    .def("add", static_cast<void (Constraints::*)(const std::string&, ConstraintComponentBasePtr)>(&Constraints::add),
          py::arg("name"), py::arg("constraint"))
    .def("add", static_cast<void (Constraints::*)(const std::string&, ImpactConstraintComponentBasePtr)>(&Constraints::add),
          py::arg("name"), py::arg("constraint"))
    .def("erase", &Constraints::exist,
          py::arg("name"))
    .def("get_constraint_component", &Constraints::getConstraintComponent,
          py::arg("name"))
    .def("get_impact_constraint_component", &Constraints::getImpactConstraintComponent,
          py::arg("name"))
    .def("clear", &Constraints::clear)
    .def("set_barrier_param", &Constraints::setBarrierParam,
          py::arg("barrier_param"))
    .def("set_fraction_to_boundary_rule", &Constraints::setFractionToBoundaryRule,
          py::arg("fraction_to_boundary_rule"))
    .def("get_barrier_param", &Constraints::getBarrierParam)
    .def("get_fraction_to_boundary_rule", &Constraints::getFractionToBoundaryRule)
    DEFINE_ROBOTOC_PYBIND11_CLASS_CLONE(Constraints);
}

} // namespace python
} // namespace robotoc