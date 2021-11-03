#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "robotoc/hybrid/sto_constraints.hpp"


namespace robotoc {
namespace python {

namespace py = pybind11;

PYBIND11_MODULE(sto_constraints, m) {
  py::class_<STOConstraints, std::shared_ptr<STOConstraints>>(m, "STOConstraints")
    .def(py::init<>())
    .def("set_barrier", &STOConstraints::setBarrier)
    .def("set_fraction_to_boundary_rule", &STOConstraints::setFractionToBoundaryRule)
    .def("set_minimum_dwell_times", static_cast<void (STOConstraints::*)(const double)>(&STOConstraints::setMinimumDwellTimes))
    .def("set_minimum_dwell_times", static_cast<void (STOConstraints::*)(const std::vector<double>&)>(&STOConstraints::setMinimumDwellTimes));
}

} // namespace python
} // namespace robotoc