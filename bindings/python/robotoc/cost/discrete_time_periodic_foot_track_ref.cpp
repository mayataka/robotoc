#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include <pybind11/numpy.h>

#include "robotoc/cost/discrete_time_periodic_foot_track_ref.hpp"


namespace robotoc {
namespace python {

namespace py = pybind11;

PYBIND11_MODULE(discrete_time_periodic_foot_track_ref, m) {
  py::class_<DiscreteTimePeriodicFootTrackRef, TimeVaryingTaskSpace3DRefBase,
             std::shared_ptr<DiscreteTimePeriodicFootTrackRef>>(m, "DiscreteTimePeriodicFootTrackRef")
    .def(py::init<const Eigen::Vector3d&, const double, const double, 
                  const int, const int, const int, const int, const bool>(),
          py::arg("x3d0"), py::arg("step_length"), py::arg("step_height"),
          py::arg("start_phase"), py::arg("end_phase"),
          py::arg("active_phases"), py::arg("inactive_phases"),
          py::arg("is_first_move_half"))
    .def("update_x3d_ref", &DiscreteTimePeriodicFootTrackRef::update_x3d_ref,
          py::arg("grid_info"), py::arg("x3d_ref"))
    .def("is_active", &DiscreteTimePeriodicFootTrackRef::isActive,
          py::arg("grid_info"));
}

} // namespace python
} // namespace robotoc