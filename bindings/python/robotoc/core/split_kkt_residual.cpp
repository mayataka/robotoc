#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/eigen.h>
#include <pybind11/numpy.h>

#include "robotoc/core/split_kkt_residual.hpp"
#include "robotoc/utils/pybind11_macros.hpp"


namespace robotoc {
namespace python {

namespace py = pybind11;

PYBIND11_MODULE(split_kkt_residual, m) {
  py::class_<SplitKKTResidual>(m, "SplitKKTResidual")
    .def(py::init<const Robot&>())
    .def(py::init<>())
    .def("set_contact_dimension", &SplitKKTResidual::setContactDimension,
          py::arg("dimf"))
    .def("set_switching_constraint_dimension", &SplitKKTResidual::setSwitchingConstraintDimension,
          py::arg("dims"))
    .def_readwrite("Fx", &SplitKKTResidual::Fx)
    .def_property("Fq", static_cast<const Eigen::VectorBlock<const Eigen::VectorXd> (SplitKKTResidual::*)() const>(&SplitKKTResidual::Fq),
                        static_cast<Eigen::VectorBlock<Eigen::VectorXd> (SplitKKTResidual::*)()>(&SplitKKTResidual::Fq))
    .def_property("Fv", static_cast<const Eigen::VectorBlock<const Eigen::VectorXd> (SplitKKTResidual::*)() const>(&SplitKKTResidual::Fv),
                        static_cast<Eigen::VectorBlock<Eigen::VectorXd> (SplitKKTResidual::*)()>(&SplitKKTResidual::Fv))
    .def_property("P", static_cast<const Eigen::VectorBlock<const Eigen::VectorXd> (SplitKKTResidual::*)() const>(&SplitKKTResidual::P),
                       static_cast<Eigen::VectorBlock<Eigen::VectorXd> (SplitKKTResidual::*)()>(&SplitKKTResidual::P))
    .def_readwrite("lx", &SplitKKTResidual::lx)
    .def_property("lq", static_cast<const Eigen::VectorBlock<const Eigen::VectorXd> (SplitKKTResidual::*)() const>(&SplitKKTResidual::lq),
                        static_cast<Eigen::VectorBlock<Eigen::VectorXd> (SplitKKTResidual::*)()>(&SplitKKTResidual::lq))
    .def_property("lv", static_cast<const Eigen::VectorBlock<const Eigen::VectorXd> (SplitKKTResidual::*)() const>(&SplitKKTResidual::lv),
                        static_cast<Eigen::VectorBlock<Eigen::VectorXd> (SplitKKTResidual::*)()>(&SplitKKTResidual::lv))
    .def_readwrite("la", &SplitKKTResidual::la)
    .def_readwrite("lu", &SplitKKTResidual::lu)
    .def_readwrite("h", &SplitKKTResidual::h)
    .def_property("lf", static_cast<const Eigen::VectorBlock<const Eigen::VectorXd> (SplitKKTResidual::*)() const>(&SplitKKTResidual::lf),
                        static_cast<Eigen::VectorBlock<Eigen::VectorXd> (SplitKKTResidual::*)()>(&SplitKKTResidual::lf))
    .def("dimf", &SplitKKTResidual::dimf)
    .def("dims", &SplitKKTResidual::dims)
    DEFINE_ROBOTOC_PYBIND11_CLASS_CLONE(SplitKKTResidual)
    DEFINE_ROBOTOC_PYBIND11_CLASS_PRINT(SplitKKTResidual);
}

} // namespace python
} // namespace robotoc