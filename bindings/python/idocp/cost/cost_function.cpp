#include <pybind11/pybind11.h>

#include "idocp/cost/cost_function.hpp"


namespace idocp {
namespace python {

namespace py = pybind11;

PYBIND11_MODULE(cost_function, m) {
  py::class_<CostFunction, std::shared_ptr<CostFunction>>(m, "CostFunction")
    .def(py::init<>())
    .def("push_back", &CostFunction::push_back)
    .def("clear", &CostFunction::clear);

  m.def("create_cost_function", []() {
    return std::make_shared<CostFunction>();
  });
}

} // namespace python
} // namespace idocp