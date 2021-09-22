#include <pybind11/pybind11.h>

#include "idocp/hybrid/switching_time_cost_function.hpp"


namespace idocp {
namespace python {

namespace py = pybind11;

PYBIND11_MODULE(switching_time_cost_function, m) {
  py::class_<SwitchingTimeCostFunction, std::shared_ptr<SwitchingTimeCostFunction>>(m, "SwitchingTimeCostFunction")
    .def(py::init<>())
    .def("push_back", &SwitchingTimeCostFunction::push_back)
    .def("clear", &SwitchingTimeCostFunction::clear);
}

} // namespace python
} // namespace idocp