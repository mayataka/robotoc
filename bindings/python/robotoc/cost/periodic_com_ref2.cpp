#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include <pybind11/numpy.h>

#include "robotoc/cost/periodic_com_ref2.hpp"


namespace robotoc {
namespace python {

namespace py = pybind11;

PYBIND11_MODULE(periodic_com_ref2, m) {
  py::class_<PeriodicCoMRef2, TimeVaryingCoMRefBase,
             std::shared_ptr<PeriodicCoMRef2>>(m, "PeriodicCoMRef2")
    .def(py::init<const Eigen::Vector3d&, const Eigen::Vector3d&,
                  const double, const double, const double, const bool>())
    .def("update_com_ref", &PeriodicCoMRef2::update_CoM_ref)
    .def("is_active", &PeriodicCoMRef2::isActive);
}

} // namespace python
} // namespace robotoc