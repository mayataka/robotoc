#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/eigen.h>
#include <pybind11/numpy.h>

#include "robotoc/core/split_kkt_matrix.hpp"


namespace robotoc {
namespace python {

namespace py = pybind11;

PYBIND11_MODULE(split_kkt_matrix, m) {
  py::class_<SplitKKTMatrix>(m, "SplitKKTMatrix")
    .def(py::init<const Robot&>())
    .def(py::init<>())
    // .def("set_contact_status", &SplitKKTMatrix::setContactStatus,
    //       py::arg("contact_status"))
    .def("is_dimension_consistent", &SplitKKTMatrix::isDimensionConsistent)
    .def_readwrite("Fxx", &SplitKKTMatrix::Fxx)
    .def_readwrite("Fvu", &SplitKKTMatrix::Fvu)
    .def_readwrite("Qxx", &SplitKKTMatrix::Qxx)
    .def_readwrite("Qxu", &SplitKKTMatrix::Qxu)
    .def_readwrite("Quu", &SplitKKTMatrix::Quu)
    .def_readwrite("Qaa", &SplitKKTMatrix::Qaa)
    .def_readwrite("Fqq_prev", &SplitKKTMatrix::Fqq_prev)
    .def_property("Fqq", static_cast<const Eigen::Block<const Eigen::MatrixXd> (SplitKKTMatrix::*)() const>(&SplitKKTMatrix::Fqq),
                         static_cast<Eigen::Block<Eigen::MatrixXd> (SplitKKTMatrix::*)()>(&SplitKKTMatrix::Fqq))
    .def_property("Fqv", static_cast<const Eigen::Block<const Eigen::MatrixXd> (SplitKKTMatrix::*)() const>(&SplitKKTMatrix::Fqv),
                         static_cast<Eigen::Block<Eigen::MatrixXd> (SplitKKTMatrix::*)()>(&SplitKKTMatrix::Fqv))
    .def_property("Fvq", static_cast<const Eigen::Block<const Eigen::MatrixXd> (SplitKKTMatrix::*)() const>(&SplitKKTMatrix::Fvq),
                         static_cast<Eigen::Block<Eigen::MatrixXd> (SplitKKTMatrix::*)()>(&SplitKKTMatrix::Fvq))
    .def_property("Fvv", static_cast<const Eigen::Block<const Eigen::MatrixXd> (SplitKKTMatrix::*)() const>(&SplitKKTMatrix::Fvv),
                         static_cast<Eigen::Block<Eigen::MatrixXd> (SplitKKTMatrix::*)()>(&SplitKKTMatrix::Fvv))
    .def_property("Qqq", static_cast<const Eigen::Block<const Eigen::MatrixXd> (SplitKKTMatrix::*)() const>(&SplitKKTMatrix::Qqq),
                         static_cast<Eigen::Block<Eigen::MatrixXd> (SplitKKTMatrix::*)()>(&SplitKKTMatrix::Qqq))
    .def_property("Qqv", static_cast<const Eigen::Block<const Eigen::MatrixXd> (SplitKKTMatrix::*)() const>(&SplitKKTMatrix::Qqv),
                         static_cast<Eigen::Block<Eigen::MatrixXd> (SplitKKTMatrix::*)()>(&SplitKKTMatrix::Qqv))
    .def_property("Qvq", static_cast<const Eigen::Block<const Eigen::MatrixXd> (SplitKKTMatrix::*)() const>(&SplitKKTMatrix::Qvq),
                         static_cast<Eigen::Block<Eigen::MatrixXd> (SplitKKTMatrix::*)()>(&SplitKKTMatrix::Qvq))
    .def_property("Qvv", static_cast<const Eigen::Block<const Eigen::MatrixXd> (SplitKKTMatrix::*)() const>(&SplitKKTMatrix::Qvv),
                         static_cast<Eigen::Block<Eigen::MatrixXd> (SplitKKTMatrix::*)()>(&SplitKKTMatrix::Qvv))
    .def_property("Qqu", static_cast<const Eigen::Block<const Eigen::MatrixXd> (SplitKKTMatrix::*)() const>(&SplitKKTMatrix::Qqu),
                         static_cast<Eigen::Block<Eigen::MatrixXd> (SplitKKTMatrix::*)()>(&SplitKKTMatrix::Qqu))
    .def_property("Qvu", static_cast<const Eigen::Block<const Eigen::MatrixXd> (SplitKKTMatrix::*)() const>(&SplitKKTMatrix::Qvu),
                         static_cast<Eigen::Block<Eigen::MatrixXd> (SplitKKTMatrix::*)()>(&SplitKKTMatrix::Qvu))
    .def_property("Qff", static_cast<const Eigen::Block<const Eigen::MatrixXd> (SplitKKTMatrix::*)() const>(&SplitKKTMatrix::Qff),
                         static_cast<Eigen::Block<Eigen::MatrixXd> (SplitKKTMatrix::*)()>(&SplitKKTMatrix::Qff))
    .def_property("Qqf", static_cast<const Eigen::Block<const Eigen::MatrixXd> (SplitKKTMatrix::*)() const>(&SplitKKTMatrix::Qqf),
                         static_cast<Eigen::Block<Eigen::MatrixXd> (SplitKKTMatrix::*)()>(&SplitKKTMatrix::Qqf))
    .def("__str__", [](const SplitKKTMatrix& self) {
        std::stringstream ss;
        ss << self;
        return ss.str();
      });
}

} // namespace python
} // namespace robotoc