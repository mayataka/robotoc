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
                  const double, const double, const double, const bool>(),
          py::arg("com_ref0"), py::arg("vcom_ref"), py::arg("t0"),
          py::arg("period_active"), py::arg("period_inactive"), 
          py::arg("is_first_move_half"))
    .def("update_com_ref", &PeriodicCoMRef2::update_com_ref,
          py::arg("t"), py::arg("com_ref"))
    .def("is_active", &PeriodicCoMRef2::isActive,
          py::arg("t"));
}

} // namespace python
} // namespace robotoc