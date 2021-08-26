#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include <pybind11/numpy.h>

#include "idocp/cost/periodic_foot_track_ref2.hpp"


namespace idocp {
namespace python {

namespace py = pybind11;

PYBIND11_MODULE(periodic_foot_track_ref2, m) {
  py::class_<PeriodicFootTrackRef2, TimeVaryingTaskSpace3DRefBase,
             std::shared_ptr<PeriodicFootTrackRef2>>(m, "PeriodicFootTrackRef2")
    .def(py::init<const Eigen::Vector3d&, const double, const double, 
                  const double, const double, const double, const bool>())
    .def("update_q_3d_ref", &PeriodicFootTrackRef2::update_q_3d_ref)
    .def("is_active", &PeriodicFootTrackRef2::isActive);
}

} // namespace python
} // namespace idocp