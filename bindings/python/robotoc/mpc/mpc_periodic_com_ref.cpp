#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>

#include "robotoc/mpc/mpc_periodic_com_ref.hpp"


namespace robotoc {
namespace python {

namespace py = pybind11;

PYBIND11_MODULE(mpc_periodic_com_ref, m) {
  py::class_<MPCPeriodicCoMRef, TimeVaryingCoMRefBase,
             std::shared_ptr<MPCPeriodicCoMRef>>(m, "MPCPeriodicCoMRef")
    .def(py::init<const double, const double, const double>(),
          py::arg("swing_start_time"), py::arg("period_active"), 
          py::arg("period_inactive"))
    .def("set_period", &MPCPeriodicCoMRef::setPeriod,
          py::arg("swing_start_time"), py::arg("period_active"), 
          py::arg("period_inactive"))
    .def("set_com_ref", &MPCPeriodicCoMRef::setCoMRef,
          py::arg("contact_sequence"), py::arg("foot_step_planner"))
    .def("update_com_ref", &MPCPeriodicCoMRef::update_com_ref,
          py::arg("grid_info"), py::arg("com_ref"))
    .def("is_active", &MPCPeriodicCoMRef::isActive,
          py::arg("grid_info"));
}

} // namespace python
} // namespace robotoc