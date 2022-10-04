#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/eigen.h>
#include <pybind11/numpy.h>

#include "robotoc/robot/contact_model_info.hpp"
#include "robotoc/utils/pybind11_macros.hpp"

#include <iostream>


namespace robotoc {
namespace python {

namespace py = pybind11;

PYBIND11_MODULE(contact_model_info, m) {
  py::class_<ContactModelInfo>(m, "ContactModelInfo")
    .def(py::init<const std::string&, const double>(), 
          py::arg("frame"), py::arg("baumgarte_time_step"))
    .def(py::init<const std::string&, const double, const double>(), 
          py::arg("frame"), py::arg("baumgarte_position_gain"), 
          py::arg("baumgarte_velocity_gain"))
    .def(py::init<>())
    .def_readwrite("frame", &ContactModelInfo::frame)
    .def_readwrite("baumgarte_position_gain", &ContactModelInfo::baumgarte_position_gain)
    .def_readwrite("baumgarte_velocity_gain", &ContactModelInfo::baumgarte_velocity_gain) 
     DEFINE_ROBOTOC_PYBIND11_CLASS_CLONE(ContactModelInfo);
}

} // namespace python
} // namespace robotoc