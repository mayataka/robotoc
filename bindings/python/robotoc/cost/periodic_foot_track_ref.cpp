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
                  const double, const double, const double, const bool>())
    .def("update_q_3d_ref", &PeriodicFootTrackRef::update_q_3d_ref)
    .def("is_active", &PeriodicFootTrackRef::isActive);
}

} // namespace python
} // namespace robotoc