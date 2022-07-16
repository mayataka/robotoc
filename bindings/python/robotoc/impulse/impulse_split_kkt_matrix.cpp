#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/eigen.h>
#include <pybind11/numpy.h>

#include "robotoc/impulse/impulse_split_kkt_matrix.hpp"


namespace robotoc {
namespace python {

namespace py = pybind11;

PYBIND11_MODULE(impulse_split_kkt_matrix, m) {
  py::class_<ImpulseSplitKKTMatrix>(m, "ImpulseSplitKKTMatrix")
    .def(py::init<const Robot&>())
    .def(py::init<>())
    .def("clone", [](const ImpulseSplitKKTMatrix& self) {
       auto other = self;
       return other;
     })
    .def("set_impulse_status", &ImpulseSplitKKTMatrix::setImpulseStatus,
          py::arg("impulse_status"))
    .def("is_dimension_consistent", &ImpulseSplitKKTMatrix::isDimensionConsistent)
    .def_readwrite("Fxx", &ImpulseSplitKKTMatrix::Fxx)
    .def_readwrite("Qxx", &ImpulseSplitKKTMatrix::Qxx)
    .def_readwrite("Qdvdvd", &ImpulseSplitKKTMatrix::Qdvdv)
    .def_readwrite("Fqq_prev", &ImpulseSplitKKTMatrix::Fqq_prev)
    .def_property("Fqq", static_cast<const Eigen::Block<const Eigen::MatrixXd> (ImpulseSplitKKTMatrix::*)() const>(&ImpulseSplitKKTMatrix::Fqq),
                         static_cast<Eigen::Block<Eigen::MatrixXd> (ImpulseSplitKKTMatrix::*)()>(&ImpulseSplitKKTMatrix::Fqq))
    .def_property("Fqv", static_cast<const Eigen::Block<const Eigen::MatrixXd> (ImpulseSplitKKTMatrix::*)() const>(&ImpulseSplitKKTMatrix::Fqv),
                         static_cast<Eigen::Block<Eigen::MatrixXd> (ImpulseSplitKKTMatrix::*)()>(&ImpulseSplitKKTMatrix::Fqv))
    .def_property("Fvq", static_cast<const Eigen::Block<const Eigen::MatrixXd> (ImpulseSplitKKTMatrix::*)() const>(&ImpulseSplitKKTMatrix::Fvq),
                         static_cast<Eigen::Block<Eigen::MatrixXd> (ImpulseSplitKKTMatrix::*)()>(&ImpulseSplitKKTMatrix::Fvq))
    .def_property("Fvv", static_cast<const Eigen::Block<const Eigen::MatrixXd> (ImpulseSplitKKTMatrix::*)() const>(&ImpulseSplitKKTMatrix::Fvv),
                         static_cast<Eigen::Block<Eigen::MatrixXd> (ImpulseSplitKKTMatrix::*)()>(&ImpulseSplitKKTMatrix::Fvv))
    .def_property("Qqq", static_cast<const Eigen::Block<const Eigen::MatrixXd> (ImpulseSplitKKTMatrix::*)() const>(&ImpulseSplitKKTMatrix::Qqq),
                         static_cast<Eigen::Block<Eigen::MatrixXd> (ImpulseSplitKKTMatrix::*)()>(&ImpulseSplitKKTMatrix::Qqq))
    .def_property("Qqv", static_cast<const Eigen::Block<const Eigen::MatrixXd> (ImpulseSplitKKTMatrix::*)() const>(&ImpulseSplitKKTMatrix::Qqv),
                         static_cast<Eigen::Block<Eigen::MatrixXd> (ImpulseSplitKKTMatrix::*)()>(&ImpulseSplitKKTMatrix::Qqv))
    .def_property("Qvq", static_cast<const Eigen::Block<const Eigen::MatrixXd> (ImpulseSplitKKTMatrix::*)() const>(&ImpulseSplitKKTMatrix::Qvq),
                         static_cast<Eigen::Block<Eigen::MatrixXd> (ImpulseSplitKKTMatrix::*)()>(&ImpulseSplitKKTMatrix::Qvq))
    .def_property("Qvv", static_cast<const Eigen::Block<const Eigen::MatrixXd> (ImpulseSplitKKTMatrix::*)() const>(&ImpulseSplitKKTMatrix::Qvv),
                         static_cast<Eigen::Block<Eigen::MatrixXd> (ImpulseSplitKKTMatrix::*)()>(&ImpulseSplitKKTMatrix::Qvv))
    .def_property("Qff", static_cast<const Eigen::Block<const Eigen::MatrixXd> (ImpulseSplitKKTMatrix::*)() const>(&ImpulseSplitKKTMatrix::Qff),
                         static_cast<Eigen::Block<Eigen::MatrixXd> (ImpulseSplitKKTMatrix::*)()>(&ImpulseSplitKKTMatrix::Qff))
    .def_property("Qqf", static_cast<const Eigen::Block<const Eigen::MatrixXd> (ImpulseSplitKKTMatrix::*)() const>(&ImpulseSplitKKTMatrix::Qqf),
                         static_cast<Eigen::Block<Eigen::MatrixXd> (ImpulseSplitKKTMatrix::*)()>(&ImpulseSplitKKTMatrix::Qqf))
    .def("__str__", [](const ImpulseSplitKKTMatrix& self) {
        std::stringstream ss;
        ss << self;
        return ss.str();
      });
}

} // namespace python
} // namespace robotoc