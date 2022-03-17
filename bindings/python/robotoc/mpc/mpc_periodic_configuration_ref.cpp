#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>

#include "robotoc/mpc/mpc_periodic_configuration_ref.hpp"


namespace robotoc {
namespace python {

namespace py = pybind11;

PYBIND11_MODULE(mpc_periodic_configuration_ref, m) {
  py::class_<MPCPeriodicConfigurationRef, TimeVaryingConfigurationRefBase,
             std::shared_ptr<MPCPeriodicConfigurationRef>>(m, "MPCPeriodicConfigurationRef")
    .def(py::init<const Eigen::VectorXd&, const double, const double, const double>(),
          py::arg("q"), py::arg("swing_start_time"), py::arg("period_active"), 
          py::arg("period_inactive"))
    .def("set_period", &MPCPeriodicConfigurationRef::setPeriod,
          py::arg("swing_start_time"), py::arg("period_active"), 
          py::arg("period_inactive"))
    .def("set_configuration_ref", &MPCPeriodicConfigurationRef::setConfigurationRef,
          py::arg("contact_sequence"), py::arg("foot_step_planner"))
    .def("update_q_ref", &MPCPeriodicConfigurationRef::update_q_ref,
          py::arg("robot"), py::arg("grid_info"), py::arg("com_ref"))
    .def("is_active", &MPCPeriodicConfigurationRef::isActive,
          py::arg("grid_info"));
}

} // namespace python
} // namespace robotoc