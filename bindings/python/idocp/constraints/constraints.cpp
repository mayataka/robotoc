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
    .def("append", static_cast<void (Constraints::*)(const ConstraintComponentBasePtr&)>(&Constraints::push_back))
    .def("append", static_cast<void (Constraints::*)(const ImpulseConstraintComponentBasePtr&)>(&Constraints::push_back))
    .def("clear", &Constraints::clear);

  m.def("create_constraints", []() {
    return std::make_shared<Constraints>();
  });
}

} // namespace python
} // namespace idocp