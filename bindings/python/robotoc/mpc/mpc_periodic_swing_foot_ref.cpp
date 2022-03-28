#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>

#include "robotoc/mpc/mpc_periodic_swing_foot_ref.hpp"


namespace robotoc {
namespace python {

namespace py = pybind11;

PYBIND11_MODULE(mpc_periodic_swing_foot_ref, m) {
  py::class_<MPCPeriodicSwingFootRef, TimeVaryingTaskSpace3DRefBase,
             std::shared_ptr<MPCPeriodicSwingFootRef>>(m, "MPCPeriodicSwingFootRef")
    .def(py::init<const int, const double, const double, const double, const double>(),
          py::arg("contact_index"), py::arg("swing_height"), 
          py::arg("swing_start_time"), py::arg("period_active"), py::arg("period_inactive"))
    .def("set_period", &MPCPeriodicSwingFootRef::setPeriod,
          py::arg("swing_start_time"), py::arg("period_active"), 
          py::arg("period_inactive"))
    .def("set_swing_foot_ref", &MPCPeriodicSwingFootRef::setSwingFootRef,
          py::arg("contact_sequence"), py::arg("foot_step_planner"))
    .def("update_x3d_ref", &MPCPeriodicSwingFootRef::update_x3d_ref,
          py::arg("grid_info"), py::arg("x3d_ref"))
    .def("is_active", &MPCPeriodicSwingFootRef::isActive,
          py::arg("grid_info"));
}

} // namespace python
} // namespace robotoc