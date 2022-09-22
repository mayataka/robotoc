#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include <pybind11/numpy.h>

#include "robotoc/sto/sto_cost_function.hpp"


namespace robotoc {
namespace python {

namespace py = pybind11;

PYBIND11_MODULE(sto_cost_function, m) {
  py::class_<STOCostFunction, std::shared_ptr<STOCostFunction>>(m, "STOCostFunction")
    .def(py::init<>())
    .def("push_back", &STOCostFunction::push_back)
    .def("clear", &STOCostFunction::clear);
}

} // namespace python
} // namespace robotoc