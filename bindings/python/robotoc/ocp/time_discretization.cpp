#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/eigen.h>
#include <pybind11/numpy.h>

#include "robotoc/ocp/discretization_method.hpp"
#include "robotoc/ocp/time_discretization.hpp"

#include <iostream>


namespace robotoc {
namespace python {

namespace py = pybind11;

PYBIND11_MODULE(time_discretization, m) {
  py::enum_<DiscretizationMethod>(m, "DiscretizationMethod", py::arithmetic())
    .value("GridBased",  DiscretizationMethod::GridBased)
    .value("PhaseBased", DiscretizationMethod::PhaseBased)
    .export_values();

  py::class_<TimeDiscretization>(m, "TimeDiscretization")
    .def(py::init<const double, const int, const int>(), 
          py::arg("T"), py::arg("N"), py::arg("reserved_num_discrete_events")=0)
    .def("set_discretization_method", &TimeDiscretization::setDiscretizationMethod,
          py::arg("discretization_method"))
    .def("N", &TimeDiscretization::N)
    .def("discretizeGrid", &TimeDiscretization::discretizeGrid,
          py::arg("contact_sequence"), py::arg("t")) 
    .def("discretizePhase", &TimeDiscretization::discretizePhase,
          py::arg("contact_sequence"), py::arg("t")) 
    .def("getGrid", &TimeDiscretization::getGrid) 
    .def("N_grids", &TimeDiscretization::N_grids) 
    .def("__str__", [](const TimeDiscretization& self) {
        std::stringstream ss;
        ss << self;
        return ss.str();
      });
}

} // namespace python
} // namespace robotoc