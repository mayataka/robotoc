#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include <pybind11/numpy.h>

#include "robotoc/cost/periodic_swing_foot_ref.hpp"


namespace robotoc {
namespace python {

namespace py = pybind11;

PYBIND11_MODULE(periodic_swing_foot_ref, m) {
  py::class_<PeriodicSwingFootRef, TaskSpace3DRefBase,
             std::shared_ptr<PeriodicSwingFootRef>>(m, "PeriodicSwingFootRef")
    .def(py::init<const Eigen::Vector3d&, const Eigen::Vector3d&, const double, 
                  const double, const double, const double, const bool>(),
          py::arg("x3d0"), py::arg("step_length"), py::arg("step_height"),
          py::arg("t0"), py::arg("period_swing"), py::arg("period_stance"),
          py::arg("is_first_step_half"))
    .def("clone", [](const PeriodicSwingFootRef& self) {
       auto other = self;
       return other;
     })
    .def("set_foot_track_ref", &PeriodicSwingFootRef::setFootTrackRef,
          py::arg("x3d0"), py::arg("step_length"), py::arg("step_height"),
          py::arg("t0"), py::arg("period_swing"), py::arg("period_stance"),
          py::arg("is_first_step_half")=false)
    .def("updateRef", &PeriodicSwingFootRef::updateRef,
          py::arg("grid_info"), py::arg("x3d_ref"))
    .def("is_active", &PeriodicSwingFootRef::isActive,
          py::arg("grid_info"));
}

} // namespace python
} // namespace robotoc