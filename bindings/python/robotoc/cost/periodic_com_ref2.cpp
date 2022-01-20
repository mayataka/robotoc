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
                  const int, const int, const int, const int, const bool>(),
          py::arg("com_ref0"), py::arg("com_step"), 
          py::arg("start_phase"), py::arg("end_phase"),
          py::arg("active_phases"), py::arg("inactive_phases"),
          py::arg("is_first_move_half"))
    .def("update_com_ref", &PeriodicCoMRef2::update_com_ref,
          py::arg("grid_info"), py::arg("com_ref"))
    .def("is_active", &PeriodicCoMRef2::isActive,
          py::arg("grid_info"));
}

} // namespace python
} // namespace robotoc