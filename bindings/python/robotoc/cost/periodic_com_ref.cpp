#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include <pybind11/numpy.h>

#include "robotoc/cost/periodic_com_ref.hpp"


namespace robotoc {
namespace python {

namespace py = pybind11;

PYBIND11_MODULE(periodic_com_ref, m) {
  py::class_<PeriodicCoMRef, TimeVaryingCoMRefBase,
             std::shared_ptr<PeriodicCoMRef>>(m, "PeriodicCoMRef")
    .def(py::init<const Eigen::Vector3d&, const Eigen::Vector3d&,
                  const double, const double, const double, const bool>())
    .def("update_com_ref", &PeriodicCoMRef::update_CoM_ref)
    .def("is_active", &PeriodicCoMRef::isActive);
}

} // namespace python
} // namespace robotoc