#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include <pybind11/numpy.h>

#include "robotoc/cost/cost_function_data.hpp"
#include "robotoc/utils/pybind11_macros.hpp"


namespace robotoc {
namespace python {

namespace py = pybind11;

PYBIND11_MODULE(cost_function_data, m) {
  py::class_<CostFunctionData>(m, "CostFunctionData")
    .def(py::init<Robot>(),
          py::arg("robot"))
    .def(py::init<>())
    .def_readwrite("qdiff", &CostFunctionData::qdiff)
    .def_readwrite("q_ref", &CostFunctionData::q_ref)
    .def_readwrite("x3d_ref", &CostFunctionData::x3d_ref)
    .def_readwrite("diff_3d", &CostFunctionData::diff_3d)
    .def_readwrite("diff_6d", &CostFunctionData::diff_6d)
    .def_readwrite("x6d_ref", &CostFunctionData::x6d_ref)
    .def_readwrite("x6d_ref_inv", &CostFunctionData::x6d_ref_inv)
    .def_readwrite("diff_x6d", &CostFunctionData::diff_x6d)
    .def_readwrite("J_qdiff", &CostFunctionData::J_qdiff)
    .def_readwrite("J_6d", &CostFunctionData::J_6d)
    .def_readwrite("J_3d", &CostFunctionData::J_3d)
    .def_readwrite("J_66", &CostFunctionData::J_66)
    .def_readwrite("JJ_6d", &CostFunctionData::JJ_6d)
    DEFINE_ROBOTOC_PYBIND11_CLASS_CLONE(CostFunctionData);

}

} // namespace python
} // namespace robotoc