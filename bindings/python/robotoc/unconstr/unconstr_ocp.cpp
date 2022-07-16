#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "robotoc/unconstr/unconstr_ocp.hpp"


namespace robotoc {
namespace python {

namespace py = pybind11;

PYBIND11_MODULE(unconstr_ocp, m) {
  py::class_<UnconstrOCP>(m, "UnconstrOCP")
    .def(py::init<const Robot&, const std::shared_ptr<CostFunction>&,
                  const std::shared_ptr<Constraints>&, const double, const int>(),
          py::arg("robot"), py::arg("cost"), py::arg("constraints"), 
          py::arg("T"), py::arg("N"))
    .def(py::init<>())
    .def("clone", [](const UnconstrOCP& self) {
       auto other = self;
       return other;
     })
    .def("robot", &UnconstrOCP::robot)
    .def("cost", &UnconstrOCP::cost)
    .def("constraints", &UnconstrOCP::constraints)
    .def("T", &UnconstrOCP::T)
    .def("N", &UnconstrOCP::N)
    .def("__str__", [](const UnconstrOCP& self) {
        std::stringstream ss;
        ss << self;
        return ss.str();
      });
}

} // namespace python
} // namespace robotoc