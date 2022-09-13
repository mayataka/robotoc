#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/eigen.h>
#include <pybind11/numpy.h>

#include "robotoc/line_search/line_search_settings.hpp"


namespace robotoc {
namespace python {

namespace py = pybind11;

PYBIND11_MODULE(line_search_settings, m) {
  py::enum_<LineSearchMethod>(m, "LineSearchMethod", py::arithmetic())
    .value("Filter",  LineSearchMethod::Filter)
    .value("MeritBacktracking", LineSearchMethod::MeritBacktracking)
    .export_values();

  py::class_<LineSearchSettings>(m, "LineSearchSettings")
    .def(py::init<>())
    .def_readwrite("line_search_method", &LineSearchSettings::line_search_method)
    .def_readwrite("step_size_reduction_rate", &LineSearchSettings::step_size_reduction_rate)
    .def_readwrite("min_step_size", &LineSearchSettings::min_step_size)
    .def_readwrite("armijo_control_rate", &LineSearchSettings::armijo_control_rate)
    .def_readwrite("margin_rate", &LineSearchSettings::margin_rate)
    .def_readwrite("eps", &LineSearchSettings::eps)
    .def("__str__", [](const LineSearchSettings& self) {
        std::stringstream ss;
        ss << self;
        return ss.str();
      });
}

} // namespace python
} // namespace robotoc