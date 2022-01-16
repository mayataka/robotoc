#include <pybind11/pybind11.h>

#include "robotoc/hybrid/grid_info.hpp"

#include <iostream>


namespace robotoc {
namespace python {

namespace py = pybind11;

PYBIND11_MODULE(grid_info, m) {
  py::class_<GridInfo>(m, "GridInfo")
    .def(py::init<>())
    .def_readwrite("t", &GridInfo::t)
    .def_readwrite("dt", &GridInfo::dt)
    .def_readwrite("contact_phase", &GridInfo::contact_phase)
    .def_readwrite("time_stage", &GridInfo::time_stage)
    .def_readwrite("impulse_index", &GridInfo::impulse_index)
    .def_readwrite("lift_index", &GridInfo::lift_index)
    .def_readwrite("grid_count_in_phase", &GridInfo::grid_count_in_phase);
}

} // namespace python
} // namespace robotoc