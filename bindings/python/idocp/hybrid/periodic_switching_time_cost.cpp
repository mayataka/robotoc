#include <pybind11/pybind11.h>

#include "idocp/hybrid/periodic_switching_time_cost.hpp"


namespace idocp {
namespace python {

namespace py = pybind11;

PYBIND11_MODULE(periodic_switching_time_cost, m) {
  py::class_<PeriodicSwitchingTimeCost, SwitchingTimeCostFunctionComponentBase,
             std::shared_ptr<PeriodicSwitchingTimeCost>>(m, "PeriodicSwitchingTimeCost")
    .def(py::init<const double, const double>())
    .def("set_period", &PeriodicSwitchingTimeCost::set_period)
    .def("set_weight", &PeriodicSwitchingTimeCost::set_weight);
}

} // namespace python
} // namespace idocp