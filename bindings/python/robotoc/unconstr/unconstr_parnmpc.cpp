#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "robotoc/unconstr/unconstr_parnmpc.hpp"


namespace robotoc {
namespace python {

namespace py = pybind11;

PYBIND11_MODULE(unconstr_parnmpc, m) {
  py::class_<UnconstrParNMPC>(m, "UnconstrParNMPC")
    .def(py::init<const Robot&, const std::shared_ptr<CostFunction>&,
                  const std::shared_ptr<Constraints>&, const double, const int>(),
          py::arg("robot"), py::arg("cost"), py::arg("constraints"), 
          py::arg("T"), py::arg("N"))
    .def(py::init<>())
    .def("clone", [](const UnconstrParNMPC& self) {
       auto other = self;
       return other;
     })
    .def("robot", &UnconstrParNMPC::robot)
    .def("cost", &UnconstrParNMPC::cost)
    .def("constraints", &UnconstrParNMPC::constraints)
    .def("T", &UnconstrParNMPC::T)
    .def("N", &UnconstrParNMPC::N)
    .def("__str__", [](const UnconstrParNMPC& self) {
        std::stringstream ss;
        ss << self;
        return ss.str();
      });
}

} // namespace python
} // namespace robotoc