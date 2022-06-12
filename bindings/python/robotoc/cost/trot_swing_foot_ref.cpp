#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include <pybind11/numpy.h>

#include "robotoc/cost/trot_swing_foot_ref.hpp"


namespace robotoc {
namespace python {

namespace py = pybind11;

PYBIND11_MODULE(trot_swing_foot_ref, m) {
  py::class_<TrotSwingFootRef, SwingFootRefBase,
             std::shared_ptr<TrotSwingFootRef>>(m, "TrotSwingFootRef")
    .def(py::init<const int, const int, const int, const double, const double>(),
          py::arg("contact_index"), py::arg("x_ref_foot_contact_index"),
          py::arg("y_ref_foot_contact_index"),
          py::arg("step_length"), py::arg("step_height"))
    .def(py::init<const Robot&, const std::string&, const int, const int, const double, const double>(),
          py::arg("robot"), py::arg("contact_frame_name"), 
          py::arg("x_ref_foot_contact_index"), py::arg("y_ref_foot_contact_index"),
          py::arg("step_length"), py::arg("step_height"))
    .def("update_x3d_ref", &TrotSwingFootRef::update_x3d_ref,
          py::arg("contact_status"), py::arg("q_3d_ref"));
}

} // namespace python
} // namespace robotoc