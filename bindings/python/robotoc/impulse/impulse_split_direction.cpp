#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/eigen.h>
#include <pybind11/numpy.h>

#include "robotoc/impulse/impulse_split_direction.hpp"


namespace robotoc {
namespace python {

namespace py = pybind11;

PYBIND11_MODULE(impulse_split_direction, m) {
  py::class_<ImpulseSplitDirection>(m, "ImpulseSplitDirection")
    .def(py::init<const Robot&>())
    .def(py::init<>())
    .def("clone", [](const ImpulseSplitDirection& self) {
       auto other = self;
       return other;
     })
    .def("set_impulse_status", &ImpulseSplitDirection::setImpulseStatus,
          py::arg("impulse_status"))
    .def("set_impulse_status", &ImpulseSplitDirection::setImpulseStatus,
          py::arg("impulse_status"))
    .def("is_dimension_consistent", &ImpulseSplitDirection::isDimensionConsistent)
    .def("dimi", &ImpulseSplitDirection::dimi)
    .def_readwrite("dx", &ImpulseSplitDirection::dx)
    .def_readwrite("dlmdgmm", &ImpulseSplitDirection::dlmdgmm)
    .def_property("dq", static_cast<const Eigen::VectorBlock<const Eigen::VectorXd> (ImpulseSplitDirection::*)() const>(&ImpulseSplitDirection::dq),
                        static_cast<Eigen::VectorBlock<Eigen::VectorXd> (ImpulseSplitDirection::*)()>(&ImpulseSplitDirection::dq))
    .def_property("dv", static_cast<const Eigen::VectorBlock<const Eigen::VectorXd> (ImpulseSplitDirection::*)() const>(&ImpulseSplitDirection::dv),
                        static_cast<Eigen::VectorBlock<Eigen::VectorXd> (ImpulseSplitDirection::*)()>(&ImpulseSplitDirection::dv))
    .def_property("ddvf", static_cast<const Eigen::VectorBlock<const Eigen::VectorXd> (ImpulseSplitDirection::*)() const>(&ImpulseSplitDirection::ddvf),
                         static_cast<Eigen::VectorBlock<Eigen::VectorXd> (ImpulseSplitDirection::*)()>(&ImpulseSplitDirection::ddvf))
    .def_property("ddv", static_cast<const Eigen::VectorBlock<const Eigen::VectorXd> (ImpulseSplitDirection::*)() const>(&ImpulseSplitDirection::ddv),
                         static_cast<Eigen::VectorBlock<Eigen::VectorXd> (ImpulseSplitDirection::*)()>(&ImpulseSplitDirection::ddv))
    .def_property("df", static_cast<const Eigen::VectorBlock<const Eigen::VectorXd> (ImpulseSplitDirection::*)() const>(&ImpulseSplitDirection::df),
                        static_cast<Eigen::VectorBlock<Eigen::VectorXd> (ImpulseSplitDirection::*)()>(&ImpulseSplitDirection::df))
    .def_property("dlmd", static_cast<const Eigen::VectorBlock<const Eigen::VectorXd> (ImpulseSplitDirection::*)() const>(&ImpulseSplitDirection::dlmd),
                          static_cast<Eigen::VectorBlock<Eigen::VectorXd> (ImpulseSplitDirection::*)()>(&ImpulseSplitDirection::dlmd))
    .def_property("dgmm", static_cast<const Eigen::VectorBlock<const Eigen::VectorXd> (ImpulseSplitDirection::*)() const>(&ImpulseSplitDirection::dgmm),
                          static_cast<Eigen::VectorBlock<Eigen::VectorXd> (ImpulseSplitDirection::*)()>(&ImpulseSplitDirection::dgmm))
    .def_property("dbetamu", static_cast<const Eigen::VectorBlock<const Eigen::VectorXd> (ImpulseSplitDirection::*)() const>(&ImpulseSplitDirection::dbetamu),
                             static_cast<Eigen::VectorBlock<Eigen::VectorXd> (ImpulseSplitDirection::*)()>(&ImpulseSplitDirection::dbetamu))
    .def_property("dbeta", static_cast<const Eigen::VectorBlock<const Eigen::VectorXd> (ImpulseSplitDirection::*)() const>(&ImpulseSplitDirection::dbeta),
                           static_cast<Eigen::VectorBlock<Eigen::VectorXd> (ImpulseSplitDirection::*)()>(&ImpulseSplitDirection::dbeta))
    .def_property("dmu", static_cast<const Eigen::VectorBlock<const Eigen::VectorXd> (ImpulseSplitDirection::*)() const>(&ImpulseSplitDirection::dmu),
                         static_cast<Eigen::VectorBlock<Eigen::VectorXd> (ImpulseSplitDirection::*)()>(&ImpulseSplitDirection::dmu))
    .def("__str__", [](const ImpulseSplitDirection& self) {
        std::stringstream ss;
        ss << self;
        return ss.str();
      });
}

} // namespace python
} // namespace robotoc