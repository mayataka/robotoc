#include <pybind11/pybind11.h>

#include "robotoc/constraints/contact_wrench_cone.hpp"


namespace robotoc {
namespace python {

namespace py = pybind11;

PYBIND11_MODULE(contact_wrench_cone, m) {
  py::class_<ContactWrenchCone, ConstraintComponentBase, 
             std::shared_ptr<ContactWrenchCone>>(m, "ContactWrenchCone")
    .def(py::init<const Robot&, const double, const double>(),
          py::arg("robot"), py::arg("X"), py::arg("Y"))
    .def("set_rectangular", &ContactWrenchCone::setRectangular,
          py::arg("X"), py::arg("Y"));
}

} // namespace python
} // namespace robotoc