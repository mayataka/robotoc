#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>

#include "robotoc/cost/local_contact_force_cost.hpp"


namespace robotoc {
namespace python {

namespace py = pybind11;

PYBIND11_MODULE(local_contact_force_cost, m) {
  py::class_<LocalContactForceCost, CostFunctionComponentBase,
             std::shared_ptr<LocalContactForceCost>>(m, "LocalContactForceCost")
    .def(py::init<const Robot&>())
    .def("set_f_ref", &LocalContactForceCost::set_f_ref)
    .def("set_fi_ref", &LocalContactForceCost::set_fi_ref)
    .def("set_f_weight", &LocalContactForceCost::set_f_weight)
    .def("set_fi_weight", &LocalContactForceCost::set_fi_weight);
}

} // namespace python
} // namespace robotoc