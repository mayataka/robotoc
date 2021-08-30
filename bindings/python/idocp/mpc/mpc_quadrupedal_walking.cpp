#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/eigen.h>
#include <pybind11/numpy.h>

#include "idocp/mpc/mpc_quadrupedal_walking.hpp"


namespace idocp {
namespace python {

namespace py = pybind11;

PYBIND11_MODULE(mpc_quadrupedal_walking, m) {
  py::class_<MPCQuadrupedalWalking>(m, "MPCQuadrupedalWalking")
    .def(py::init<const Robot&, const std::shared_ptr<CostFunction>&,
                  const std::shared_ptr<Constraints>&, const double, const int, 
                  const int, const int>(),
         py::arg("robot"), py::arg("cost"), py::arg("constraints"),
         py::arg("T"), py::arg("N"), py::arg("max_num_steps"),
         py::arg("nthreads"))
    .def("set_gait_pattern", &MPCQuadrupedalWalking::setGaitPattern,
         py::arg("step_length"), py::arg("step_height"), py::arg("swing_time"), 
         py::arg("t0"))
    .def("init", &MPCQuadrupedalWalking::init,
          py::arg("t"), py::arg("q"), py::arg("v"), py::arg("num_iteration"))
    .def("update_solution", &MPCQuadrupedalWalking::updateSolution,
          py::arg("t"), py::arg("q"), py::arg("v"), py::arg("num_iteration"))
    .def("get_initial_control_input", &MPCQuadrupedalWalking::getInitialControlInput)
    .def("KKT_error", static_cast<double (MPCQuadrupedalWalking::*)()>(&MPCQuadrupedalWalking::KKTError))
    .def("KKT_error", static_cast<double (MPCQuadrupedalWalking::*)(const double, const Eigen::VectorXd&, const Eigen::VectorXd&)>(&MPCQuadrupedalWalking::KKTError),
          py::arg("t"), py::arg("q"), py::arg("v"))
//     .def("check_formulation", &MPCQuadrupedalWalking::checkFormulation)
    .def("show_info", &MPCQuadrupedalWalking::showInfo);
}

} // namespace python
} // namespace idocp