#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/eigen.h>
#include <pybind11/numpy.h>

#include "robotoc/ocp/split_kkt_residual.hpp"


namespace robotoc {
namespace python {

namespace py = pybind11;

PYBIND11_MODULE(split_kkt_residual, m) {
  py::class_<SplitKKTResidual>(m, "SplitKKTResidual")
    .def(py::init<const Robot&>())
    .def(py::init<>())
    .def("clone", [](const SplitKKTResidual& self) {
       auto other = self;
       return other;
     })
    .def("set_contact_status", &SplitKKTResidual::setContactStatus,
          py::arg("contact_status"))
    .def("is_dimension_consistent", &SplitKKTResidual::isDimensionConsistent)
    .def_readwrite("Fx", &SplitKKTResidual::Fx)
    .def_readwrite("lx", &SplitKKTResidual::lx)
    .def_readwrite("la", &SplitKKTResidual::la)
    .def_readwrite("lu", &SplitKKTResidual::lu)
    .def_readwrite("h", &SplitKKTResidual::h)
    .def_readwrite("cost", &SplitKKTResidual::cost)
    .def_property("Fq", static_cast<const Eigen::VectorBlock<const Eigen::VectorXd> (SplitKKTResidual::*)() const>(&SplitKKTResidual::Fq),
                        static_cast<Eigen::VectorBlock<Eigen::VectorXd> (SplitKKTResidual::*)()>(&SplitKKTResidual::Fq))
    .def_property("Fv", static_cast<const Eigen::VectorBlock<const Eigen::VectorXd> (SplitKKTResidual::*)() const>(&SplitKKTResidual::Fv),
                        static_cast<Eigen::VectorBlock<Eigen::VectorXd> (SplitKKTResidual::*)()>(&SplitKKTResidual::Fv))
    .def_property("lq", static_cast<const Eigen::VectorBlock<const Eigen::VectorXd> (SplitKKTResidual::*)() const>(&SplitKKTResidual::lq),
                        static_cast<Eigen::VectorBlock<Eigen::VectorXd> (SplitKKTResidual::*)()>(&SplitKKTResidual::lq))
    .def_property("lv", static_cast<const Eigen::VectorBlock<const Eigen::VectorXd> (SplitKKTResidual::*)() const>(&SplitKKTResidual::lv),
                        static_cast<Eigen::VectorBlock<Eigen::VectorXd> (SplitKKTResidual::*)()>(&SplitKKTResidual::lv))
    .def_property("lf", static_cast<const Eigen::VectorBlock<const Eigen::VectorXd> (SplitKKTResidual::*)() const>(&SplitKKTResidual::lf),
                        static_cast<Eigen::VectorBlock<Eigen::VectorXd> (SplitKKTResidual::*)()>(&SplitKKTResidual::lf))
    .def("__str__", [](const SplitKKTResidual& self) {
        std::stringstream ss;
        ss << self;
        return ss.str();
      });
}

} // namespace python
} // namespace robotoc