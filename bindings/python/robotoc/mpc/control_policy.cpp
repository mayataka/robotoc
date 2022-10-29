#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/eigen.h>
#include <pybind11/numpy.h>

#include "robotoc/mpc/control_policy.hpp"
#include "robotoc/utils/pybind11_macros.hpp"


namespace robotoc {
namespace python {

namespace py = pybind11;

PYBIND11_MODULE(control_policy, m) {
  py::class_<ControlPolicy>(m, "ControlPolicy")
    .def(py::init<const OCPSolver&, const double>(),
         py::arg("ocp_solver"), py::arg("t"))
    .def_readwrite("t", &ControlPolicy::t)
    .def_readwrite("tauJ", &ControlPolicy::tauJ)
    .def_readwrite("qJ", &ControlPolicy::qJ)
    .def_readwrite("dqJ", &ControlPolicy::dqJ)
    .def_readwrite("Kp", &ControlPolicy::Kp)
    .def_readwrite("Kd", &ControlPolicy::Kd)
    .def("set", &ControlPolicy::set,
         py::arg("ocp_solver"), py::arg("t"))
    DEFINE_ROBOTOC_PYBIND11_CLASS_CLONE(ControlPolicy)
    DEFINE_ROBOTOC_PYBIND11_CLASS_PRINT(ControlPolicy);
}

} // namespace python
} // namespace robotoc