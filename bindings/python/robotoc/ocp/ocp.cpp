#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "robotoc/ocp/ocp.hpp"


namespace robotoc {
namespace python {

namespace py = pybind11;

PYBIND11_MODULE(ocp, m) {
  py::class_<OCP>(m, "OCP")
    .def(py::init<const Robot&, const std::shared_ptr<CostFunction>&,
                  const std::shared_ptr<Constraints>&, 
                  const std::shared_ptr<STOCostFunction>&, 
                  const std::shared_ptr<STOConstraints>&, 
                  const std::shared_ptr<ContactSequence>&, 
                  const double, const int>(),
          py::arg("robot"), py::arg("cost"), py::arg("constraints"), 
          py::arg("sto_cost"), py::arg("sto_constraints"),  
          py::arg("contact_sequence"), py::arg("T"), py::arg("N"))
    .def(py::init<const Robot&, const std::shared_ptr<CostFunction>&,
                  const std::shared_ptr<Constraints>&, 
                  const std::shared_ptr<ContactSequence>&, 
                  const double, const int>(),
          py::arg("robot"), py::arg("cost"), py::arg("constraints"), 
          py::arg("contact_sequence"), py::arg("T"), py::arg("N"))
    .def(py::init<>())
    .def("set_discretization_method", &OCP::setDiscretizationMethod, 
          py::arg("discretization_method"))
    .def("discretize", &OCP::discretize, py::arg("t"))
    .def("mesh_refinement ", &OCP::meshRefinement, py::arg("t"))
    .def("discrete", &OCP::discrete)
    .def("robot", &OCP::robot)
    .def("cost", &OCP::cost)
    .def("constraints", &OCP::constraints)
    .def("sto_cost", &OCP::sto_cost)
    .def("sto_constraints", &OCP::sto_constraints)
    .def("contact_sequence", &OCP::contact_sequence)
    .def("T", &OCP::T)
    .def("N", &OCP::N)
    .def("reserved_num_discrete_events", &OCP::reservedNumDiscreteEvents)
    .def("is_sto_enabled", &OCP::isSTOEnabled)
    .def("__str__", [](const OCP& self) {
        std::stringstream ss;
        ss << self;
        return ss.str();
      });
}

} // namespace python
} // namespace robotoc