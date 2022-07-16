#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/eigen.h>
#include <pybind11/numpy.h>

#include "robotoc/impulse/impulse_split_kkt_residual.hpp"


namespace robotoc {
namespace python {

namespace py = pybind11;

PYBIND11_MODULE(impulse_split_kkt_residual, m) {
  py::class_<ImpulseSplitKKTResidual>(m, "ImpulseSplitKKTResidual")
    .def(py::init<const Robot&>())
    .def(py::init<>())
    .def("clone", [](const ImpulseSplitKKTResidual& self) {
       auto other = self;
       return other;
     })
    .def("set_impulse_status", &ImpulseSplitKKTResidual::setImpulseStatus,
          py::arg("impulse_status"))
    .def("is_dimension_consistent", &ImpulseSplitKKTResidual::isDimensionConsistent)
    .def_readwrite("Fx", &ImpulseSplitKKTResidual::Fx)
    .def_readwrite("lx", &ImpulseSplitKKTResidual::lx)
    .def_readwrite("ldv", &ImpulseSplitKKTResidual::ldv)
    .def_readwrite("cost", &ImpulseSplitKKTResidual::cost)
    .def_property("Fq", static_cast<const Eigen::VectorBlock<const Eigen::VectorXd> (ImpulseSplitKKTResidual::*)() const>(&ImpulseSplitKKTResidual::Fq),
                        static_cast<Eigen::VectorBlock<Eigen::VectorXd> (ImpulseSplitKKTResidual::*)()>(&ImpulseSplitKKTResidual::Fq))
    .def_property("Fv", static_cast<const Eigen::VectorBlock<const Eigen::VectorXd> (ImpulseSplitKKTResidual::*)() const>(&ImpulseSplitKKTResidual::Fv),
                        static_cast<Eigen::VectorBlock<Eigen::VectorXd> (ImpulseSplitKKTResidual::*)()>(&ImpulseSplitKKTResidual::Fv))
    .def_property("lq", static_cast<const Eigen::VectorBlock<const Eigen::VectorXd> (ImpulseSplitKKTResidual::*)() const>(&ImpulseSplitKKTResidual::lq),
                        static_cast<Eigen::VectorBlock<Eigen::VectorXd> (ImpulseSplitKKTResidual::*)()>(&ImpulseSplitKKTResidual::lq))
    .def_property("lv", static_cast<const Eigen::VectorBlock<const Eigen::VectorXd> (ImpulseSplitKKTResidual::*)() const>(&ImpulseSplitKKTResidual::lv),
                        static_cast<Eigen::VectorBlock<Eigen::VectorXd> (ImpulseSplitKKTResidual::*)()>(&ImpulseSplitKKTResidual::lv))
    .def_property("lf", static_cast<const Eigen::VectorBlock<const Eigen::VectorXd> (ImpulseSplitKKTResidual::*)() const>(&ImpulseSplitKKTResidual::lf),
                        static_cast<Eigen::VectorBlock<Eigen::VectorXd> (ImpulseSplitKKTResidual::*)()>(&ImpulseSplitKKTResidual::lf))
    .def("__str__", [](const ImpulseSplitKKTResidual& self) {
        std::stringstream ss;
        ss << self;
        return ss.str();
      });
}

} // namespace python
} // namespace robotoc