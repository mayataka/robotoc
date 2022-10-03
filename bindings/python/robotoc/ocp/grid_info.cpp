#include <pybind11/pybind11.h>

#include "robotoc/ocp/grid_info.hpp"

#include <iostream>
#include <sstream>


namespace robotoc {
namespace python {

namespace py = pybind11;

PYBIND11_MODULE(grid_info, m) {
  py::enum_<GridType>(m, "GridType", py::arithmetic())
    .value("Intermediate",  GridType::Intermediate)
    .value("Impulse",  GridType::Impulse)
    .value("Terminal",  GridType::Terminal)
    .export_values();

  py::class_<GridInfo>(m, "GridInfo")
    .def(py::init<>())
    .def_readwrite("type", &GridInfo::type)
    .def_readwrite("t0", &GridInfo::t0)
    .def_readwrite("t", &GridInfo::t)
    .def_readwrite("dt", &GridInfo::dt)
    .def_readwrite("phase", &GridInfo::phase)
    .def_readwrite("stage", &GridInfo::stage)
    .def_readwrite("impulse_index", &GridInfo::impulse_index)
    .def_readwrite("lift_index", &GridInfo::lift_index)
    .def_readwrite("stage_in_phase", &GridInfo::stage_in_phase)
    .def("__str__", [](const GridInfo& self) {
        std::stringstream ss;
        ss << self;
        return ss.str();
      });
}

} // namespace python
} // namespace robotoc