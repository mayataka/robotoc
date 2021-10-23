#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include <pybind11/numpy.h>

#include "robotoc/cost/contact_force_cost.hpp"


namespace robotoc {
namespace python {

namespace py = pybind11;

PYBIND11_MODULE(contact_force_cost, m) {
  py::class_<ContactForceCost, CostFunctionComponentBase,
             std::shared_ptr<ContactForceCost>>(m, "ContactForceCost")
    .def(py::init<const Robot&>())
    .def("set_f_ref", &ContactForceCost::set_f_ref)
    .def("set_fi_ref", &ContactForceCost::set_fi_ref)
    .def("set_f_weight", &ContactForceCost::set_f_weight)
    .def("set_fi_weight", &ContactForceCost::set_fi_weight);
}

} // namespace python
} // namespace robotoc