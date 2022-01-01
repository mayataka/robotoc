#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "robotoc/robot/contact_frames.hpp"

#include <iostream>


namespace robotoc {
namespace python {

namespace py = pybind11;

PYBIND11_MODULE(contact_frames, m) {
  py::class_<ContactFrames>(m, "ContactFrames")
    .def(py::init<const std::vector<int>&, const std::vector<int>&>(),
          py::arg("point_contact_frames"), 
          py::arg("surface_contact_frames"))
    .def(py::init<>())
    .def_readwrite("point_contact_frames", &ContactFrames::point_contact_frames)
    .def_readwrite("surface_contact_frames", &ContactFrames::surface_contact_frames);
}

} // namespace python
} // namespace robotoc