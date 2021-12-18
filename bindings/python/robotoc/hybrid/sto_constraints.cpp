#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>

#include "robotoc/hybrid/sto_constraints.hpp"


namespace robotoc {
namespace python {

namespace py = pybind11;

PYBIND11_MODULE(sto_constraints, m) {
  py::class_<STOConstraints, std::shared_ptr<STOConstraints>>(m, "STOConstraints")
    .def(py::init<const int, const double, const double, const double>(),
         py::arg("max_num_switches"), 
         py::arg("min_dt")=std::sqrt(std::numeric_limits<double>::epsilon()),
         py::arg("barrier")=1.0e-03, py::arg("fraction_to_boundary_rule")=0.995)
    .def(py::init<const int, const std::vector<double>&, const double, const double>(),
         py::arg("max_num_switches"), py::arg("min_dt"),
         py::arg("barrier")=1.0e-03, py::arg("fraction_to_boundary_rule")=0.995)
    .def("set_minimum_dwell_times", static_cast<void (STOConstraints::*)(const double)>(&STOConstraints::setMinimumDwellTimes),
          py::arg("min_dt"))
    .def("set_minimum_dwell_times", static_cast<void (STOConstraints::*)(const std::vector<double>&)>(&STOConstraints::setMinimumDwellTimes),
          py::arg("min_dt"))
    .def("minimum_dwell_times", &STOConstraints::minimumDwellTimes)
    .def("set_barrier", &STOConstraints::setBarrier,
          py::arg("barrier"))
    .def("set_fraction_to_boundary_rule", &STOConstraints::setFractionToBoundaryRule,
          py::arg("fraction_to_boundary_rule"))
    .def("barrier", &STOConstraints::barrier)
    .def("fraction_to_boundary_rule", &STOConstraints::fractionToBoundaryRule);
}

} // namespace python
} // namespace robotoc