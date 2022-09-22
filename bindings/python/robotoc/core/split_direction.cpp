#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/eigen.h>
#include <pybind11/numpy.h>

#include "robotoc/core/split_direction.hpp"


namespace robotoc {
namespace python {

namespace py = pybind11;

PYBIND11_MODULE(split_direction, m) {
  py::class_<SplitDirection>(m, "SplitDirection")
    .def(py::init<const Robot&>())
    .def(py::init<>())
//     .def("set_contact_status", &SplitDirection::setContactStatus,
//           py::arg("contact_status"))
//     .def("set_switching_constraint_dimension", &SplitDirection::setSwitchingConstraintDimension,
//           py::arg("impulse_status"))
    .def("is_dimension_consistent", &SplitDirection::isDimensionConsistent)
    .def("dimf", &SplitDirection::dimf)
    .def("dimi", &SplitDirection::dimi)
    .def_readwrite("dx", &SplitDirection::dx)
    .def_readwrite("du", &SplitDirection::du)
    .def_readwrite("dlmdgmm", &SplitDirection::dlmdgmm)
    .def_readwrite("dnu_passive", &SplitDirection::dnu_passive)
    .def_property("dq", static_cast<const Eigen::VectorBlock<const Eigen::VectorXd> (SplitDirection::*)() const>(&SplitDirection::dq),
                        static_cast<Eigen::VectorBlock<Eigen::VectorXd> (SplitDirection::*)()>(&SplitDirection::dq))
    .def_property("dv", static_cast<const Eigen::VectorBlock<const Eigen::VectorXd> (SplitDirection::*)() const>(&SplitDirection::dv),
                        static_cast<Eigen::VectorBlock<Eigen::VectorXd> (SplitDirection::*)()>(&SplitDirection::dv))
    .def_property("daf", static_cast<const Eigen::VectorBlock<const Eigen::VectorXd> (SplitDirection::*)() const>(&SplitDirection::daf),
                        static_cast<Eigen::VectorBlock<Eigen::VectorXd> (SplitDirection::*)()>(&SplitDirection::daf))
    .def_property("da", static_cast<const Eigen::VectorBlock<const Eigen::VectorXd> (SplitDirection::*)() const>(&SplitDirection::da),
                        static_cast<Eigen::VectorBlock<Eigen::VectorXd> (SplitDirection::*)()>(&SplitDirection::da))
    .def_property("df", static_cast<const Eigen::VectorBlock<const Eigen::VectorXd> (SplitDirection::*)() const>(&SplitDirection::df),
                        static_cast<Eigen::VectorBlock<Eigen::VectorXd> (SplitDirection::*)()>(&SplitDirection::df))
    .def_property("dlmd", static_cast<const Eigen::VectorBlock<const Eigen::VectorXd> (SplitDirection::*)() const>(&SplitDirection::dlmd),
                          static_cast<Eigen::VectorBlock<Eigen::VectorXd> (SplitDirection::*)()>(&SplitDirection::dlmd))
    .def_property("dgmm", static_cast<const Eigen::VectorBlock<const Eigen::VectorXd> (SplitDirection::*)() const>(&SplitDirection::dgmm),
                          static_cast<Eigen::VectorBlock<Eigen::VectorXd> (SplitDirection::*)()>(&SplitDirection::dgmm))
    .def_property("dbetamu", static_cast<const Eigen::VectorBlock<const Eigen::VectorXd> (SplitDirection::*)() const>(&SplitDirection::dbetamu),
                             static_cast<Eigen::VectorBlock<Eigen::VectorXd> (SplitDirection::*)()>(&SplitDirection::dbetamu))
    .def_property("dbeta", static_cast<const Eigen::VectorBlock<const Eigen::VectorXd> (SplitDirection::*)() const>(&SplitDirection::dbeta),
                           static_cast<Eigen::VectorBlock<Eigen::VectorXd> (SplitDirection::*)()>(&SplitDirection::dbeta))
    .def_property("dmu", static_cast<const Eigen::VectorBlock<const Eigen::VectorXd> (SplitDirection::*)() const>(&SplitDirection::dmu),
                         static_cast<Eigen::VectorBlock<Eigen::VectorXd> (SplitDirection::*)()>(&SplitDirection::dmu))
    .def_property("dxi", static_cast<const Eigen::VectorBlock<const Eigen::VectorXd> (SplitDirection::*)() const>(&SplitDirection::dxi),
                         static_cast<Eigen::VectorBlock<Eigen::VectorXd> (SplitDirection::*)()>(&SplitDirection::dxi))
    .def("__str__", [](const SplitDirection& self) {
        std::stringstream ss;
        ss << self;
        return ss.str();
      });
}

} // namespace python
} // namespace robotoc