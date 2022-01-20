#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include <pybind11/numpy.h>

#include "robotoc/cost/periodic_foot_track_ref.hpp"


namespace robotoc {
namespace python {

namespace py = pybind11;

PYBIND11_MODULE(periodic_foot_track_ref, m) {
  py::class_<PeriodicFootTrackRef, TimeVaryingTaskSpace3DRefBase,
             std::shared_ptr<PeriodicFootTrackRef>>(m, "PeriodicFootTrackRef")
    .def(py::init<const Eigen::Vector3d&, const double, const double, 
                  const double, const double, const double, const bool>(),
          py::arg("x3d0"), py::arg("step_length"), py::arg("step_height"),
          py::arg("t0"), py::arg("period_swing"), py::arg("period_stance"),
          py::arg("is_first_step_half"))
    .def("update_x3d_ref", &PeriodicFootTrackRef::update_x3d_ref,
          py::arg("grid_info"), py::arg("x3d_ref"))
    .def("is_active", &PeriodicFootTrackRef::isActive,
          py::arg("grid_info"));
}

} // namespace python
} // namespace robotoc