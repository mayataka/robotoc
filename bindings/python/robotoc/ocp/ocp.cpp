#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "robotoc/ocp/ocp.hpp"
#include "robotoc/utils/pybind11_macros.hpp"


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
                  const double, const int, const int>(),
          py::arg("robot"), py::arg("cost"), py::arg("constraints"), 
          py::arg("sto_cost"), py::arg("sto_constraints"),  
          py::arg("contact_sequence"), py::arg("T"), py::arg("N"), 
          py::arg("reserved_num_discrete_events")=0)
    .def(py::init<const Robot&, const std::shared_ptr<CostFunction>&,
                  const std::shared_ptr<Constraints>&, 
                  const std::shared_ptr<ContactSequence>&, 
                  const double, const int, const int>(),
          py::arg("robot"), py::arg("cost"), py::arg("constraints"), 
          py::arg("contact_sequence"), py::arg("T"), py::arg("N"),
          py::arg("reserved_num_discrete_events")=0)
    .def(py::init<const Robot&, const std::shared_ptr<CostFunction>&,
                  const std::shared_ptr<Constraints>&, 
                  const double, const int>(),
          py::arg("robot"), py::arg("cost"), py::arg("constraints"), 
          py::arg("T"), py::arg("N"))
    .def(py::init<>())
    .def_readwrite("robot", &OCP::robot)
    .def_readwrite("cost", &OCP::cost)
    .def_readwrite("constraints", &OCP::constraints)
    .def_readwrite("sto_cost", &OCP::sto_cost)
    .def_readwrite("sto_constraints", &OCP::sto_constraints)
    .def_readwrite("contact_sequence", &OCP::contact_sequence)
    .def_readwrite("T", &OCP::T)
    .def_readwrite("N", &OCP::N)
    .def_readwrite("reserved_num_discrete_events", &OCP::reserved_num_discrete_events)
    DEFINE_ROBOTOC_PYBIND11_CLASS_CLONE(OCP)
    DEFINE_ROBOTOC_PYBIND11_CLASS_PRINT(OCP);
}

} // namespace python
} // namespace robotoc