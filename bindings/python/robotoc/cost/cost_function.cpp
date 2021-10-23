#include <pybind11/pybind11.h>

#include "robotoc/cost/cost_function.hpp"


namespace robotoc {
namespace python {

namespace py = pybind11;

PYBIND11_MODULE(cost_function, m) {
  py::class_<CostFunction, std::shared_ptr<CostFunction>>(m, "CostFunction")
    .def(py::init<>())
    .def("push_back", &CostFunction::push_back)
    .def("clear", &CostFunction::clear);
}

} // namespace python
} // namespace robotoc