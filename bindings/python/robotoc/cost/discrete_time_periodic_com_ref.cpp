#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include <pybind11/numpy.h>

#include "robotoc/cost/discrete_time_periodic_com_ref.hpp"


namespace robotoc {
namespace python {

namespace py = pybind11;

PYBIND11_MODULE(discrete_time_periodic_com_ref, m) {
  py::class_<DiscreteTimePeriodicCoMRef, TimeVaryingCoMRefBase,
             std::shared_ptr<DiscreteTimePeriodicCoMRef>>(m, "DiscreteTimePeriodicCoMRef")
    .def(py::init<const Eigen::Vector3d&, const Eigen::Vector3d&,
                  const int, const int, const int, const int, const bool>(),
          py::arg("com_ref0"), py::arg("com_step"), 
          py::arg("start_phase"), py::arg("end_phase"),
          py::arg("active_phases"), py::arg("inactive_phases"),
          py::arg("is_first_move_half"))
    .def("update_com_ref", &DiscreteTimePeriodicCoMRef::update_com_ref,
          py::arg("grid_info"), py::arg("com_ref"))
    .def("is_active", &DiscreteTimePeriodicCoMRef::isActive,
          py::arg("grid_info"));
}

} // namespace python
} // namespace robotoc