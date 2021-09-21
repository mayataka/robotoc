#include <pybind11/pybind11.h>

#include "idocp/constraints/constraints.hpp"


namespace idocp {
namespace python {

namespace py = pybind11;

using ConstraintComponentBasePtr = std::shared_ptr<ConstraintComponentBase>;
using ImpulseConstraintComponentBasePtr = std::shared_ptr<ImpulseConstraintComponentBase>;

PYBIND11_MODULE(constraints, m) {
  py::class_<Constraints, std::shared_ptr<Constraints>>(m, "Constraints")
    .def(py::init<>())
    .def("push_back", static_cast<void (Constraints::*)(const ConstraintComponentBasePtr&)>(&Constraints::push_back))
    .def("push_back", static_cast<void (Constraints::*)(const ImpulseConstraintComponentBasePtr&)>(&Constraints::push_back))
    .def("clear", &Constraints::clear)
    .def("set_barrier", &Constraints::setBarrier)
    .def("set_fraction_to_boundary_rule", &Constraints::setFractionToBoundaryRule);
}

} // namespace python
} // namespace idocp