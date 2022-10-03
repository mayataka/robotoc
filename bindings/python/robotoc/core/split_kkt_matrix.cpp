#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/eigen.h>
#include <pybind11/numpy.h>

#include "robotoc/core/split_kkt_matrix.hpp"
#include "robotoc/utils/pybind11_macros.hpp"


namespace robotoc {
namespace python {

namespace py = pybind11;

PYBIND11_MODULE(split_kkt_matrix, m) {
  py::class_<SplitKKTMatrix>(m, "SplitKKTMatrix")
    .def(py::init<const Robot&>())
    .def(py::init<>())
    .def("set_contact_dimension", &SplitKKTMatrix::setContactDimension,
          py::arg("dimf"))
    .def("set_switching_constraint_dimension", &SplitKKTMatrix::setSwitchingConstraintDimension,
          py::arg("dims"))
    .def_readwrite("Fxx", &SplitKKTMatrix::Fxx)
    .def_property("Fqq", static_cast<const Eigen::Block<const Eigen::MatrixXd> (SplitKKTMatrix::*)() const>(&SplitKKTMatrix::Fqq),
                         static_cast<Eigen::Block<Eigen::MatrixXd> (SplitKKTMatrix::*)()>(&SplitKKTMatrix::Fqq))
    .def_property("Fqv", static_cast<const Eigen::Block<const Eigen::MatrixXd> (SplitKKTMatrix::*)() const>(&SplitKKTMatrix::Fqv),
                         static_cast<Eigen::Block<Eigen::MatrixXd> (SplitKKTMatrix::*)()>(&SplitKKTMatrix::Fqv))
    .def_property("Fvq", static_cast<const Eigen::Block<const Eigen::MatrixXd> (SplitKKTMatrix::*)() const>(&SplitKKTMatrix::Fvq),
                         static_cast<Eigen::Block<Eigen::MatrixXd> (SplitKKTMatrix::*)()>(&SplitKKTMatrix::Fvq))
    .def_property("Fvv", static_cast<const Eigen::Block<const Eigen::MatrixXd> (SplitKKTMatrix::*)() const>(&SplitKKTMatrix::Fvv),
                         static_cast<Eigen::Block<Eigen::MatrixXd> (SplitKKTMatrix::*)()>(&SplitKKTMatrix::Fvv))
    .def_readwrite("Fvu", &SplitKKTMatrix::Fvu)
    .def_readwrite("fx", &SplitKKTMatrix::fx)
    .def_property("fq", static_cast<const Eigen::VectorBlock<const Eigen::VectorXd> (SplitKKTMatrix::*)() const>(&SplitKKTMatrix::fq),
                        static_cast<Eigen::VectorBlock<Eigen::VectorXd> (SplitKKTMatrix::*)()>(&SplitKKTMatrix::fq))
    .def_property("fv", static_cast<const Eigen::VectorBlock<const Eigen::VectorXd> (SplitKKTMatrix::*)() const>(&SplitKKTMatrix::fv),
                        static_cast<Eigen::VectorBlock<Eigen::VectorXd> (SplitKKTMatrix::*)()>(&SplitKKTMatrix::fv))
    .def_property("Phix", static_cast<const Eigen::Block<const Eigen::MatrixXd> (SplitKKTMatrix::*)() const>(&SplitKKTMatrix::Phix),
                          static_cast<Eigen::Block<Eigen::MatrixXd> (SplitKKTMatrix::*)()>(&SplitKKTMatrix::Phix))
    .def_property("Phiq", static_cast<const Eigen::Block<const Eigen::MatrixXd> (SplitKKTMatrix::*)() const>(&SplitKKTMatrix::Phiq),
                          static_cast<Eigen::Block<Eigen::MatrixXd> (SplitKKTMatrix::*)()>(&SplitKKTMatrix::Phiq))
    .def_property("Phiv", static_cast<const Eigen::Block<const Eigen::MatrixXd> (SplitKKTMatrix::*)() const>(&SplitKKTMatrix::Phiv),
                          static_cast<Eigen::Block<Eigen::MatrixXd> (SplitKKTMatrix::*)()>(&SplitKKTMatrix::Phiv))
    .def_property("Phia", static_cast<const Eigen::Block<const Eigen::MatrixXd> (SplitKKTMatrix::*)() const>(&SplitKKTMatrix::Phia),
                          static_cast<Eigen::Block<Eigen::MatrixXd> (SplitKKTMatrix::*)()>(&SplitKKTMatrix::Phia))
    .def_property("Phiu", static_cast<const Eigen::Block<const Eigen::MatrixXd> (SplitKKTMatrix::*)() const>(&SplitKKTMatrix::Phiu),
                          static_cast<Eigen::Block<Eigen::MatrixXd> (SplitKKTMatrix::*)()>(&SplitKKTMatrix::Phiu))
    .def_property("Phit", static_cast<const Eigen::VectorBlock<const Eigen::VectorXd> (SplitKKTMatrix::*)() const>(&SplitKKTMatrix::Phit),
                          static_cast<Eigen::VectorBlock<Eigen::VectorXd> (SplitKKTMatrix::*)()>(&SplitKKTMatrix::Phit))
    .def_readwrite("Qxx", &SplitKKTMatrix::Qxx)
    .def_property("Qqq", static_cast<const Eigen::Block<const Eigen::MatrixXd> (SplitKKTMatrix::*)() const>(&SplitKKTMatrix::Qqq),
                         static_cast<Eigen::Block<Eigen::MatrixXd> (SplitKKTMatrix::*)()>(&SplitKKTMatrix::Qqq))
    .def_property("Qqv", static_cast<const Eigen::Block<const Eigen::MatrixXd> (SplitKKTMatrix::*)() const>(&SplitKKTMatrix::Qqv),
                         static_cast<Eigen::Block<Eigen::MatrixXd> (SplitKKTMatrix::*)()>(&SplitKKTMatrix::Qqv))
    .def_property("Qvq", static_cast<const Eigen::Block<const Eigen::MatrixXd> (SplitKKTMatrix::*)() const>(&SplitKKTMatrix::Qvq),
                         static_cast<Eigen::Block<Eigen::MatrixXd> (SplitKKTMatrix::*)()>(&SplitKKTMatrix::Qvq))
    .def_property("Qvv", static_cast<const Eigen::Block<const Eigen::MatrixXd> (SplitKKTMatrix::*)() const>(&SplitKKTMatrix::Qvv),
                         static_cast<Eigen::Block<Eigen::MatrixXd> (SplitKKTMatrix::*)()>(&SplitKKTMatrix::Qvv))
    .def_readwrite("Qaa", &SplitKKTMatrix::Qaa)
    .def_readwrite("Qdvdv", &SplitKKTMatrix::Qdvdv)
    .def_readwrite("Qxu", &SplitKKTMatrix::Qxu)
    .def_property("Qqu", static_cast<const Eigen::Block<const Eigen::MatrixXd> (SplitKKTMatrix::*)() const>(&SplitKKTMatrix::Qqu),
                         static_cast<Eigen::Block<Eigen::MatrixXd> (SplitKKTMatrix::*)()>(&SplitKKTMatrix::Qqu))
    .def_property("Qvu", static_cast<const Eigen::Block<const Eigen::MatrixXd> (SplitKKTMatrix::*)() const>(&SplitKKTMatrix::Qvu),
                         static_cast<Eigen::Block<Eigen::MatrixXd> (SplitKKTMatrix::*)()>(&SplitKKTMatrix::Qvu))
    .def_readwrite("Quu", &SplitKKTMatrix::Quu)
    .def_readwrite("hx", &SplitKKTMatrix::hx)
    .def_property("hq", static_cast<const Eigen::VectorBlock<const Eigen::VectorXd> (SplitKKTMatrix::*)() const>(&SplitKKTMatrix::hq),
                        static_cast<Eigen::VectorBlock<Eigen::VectorXd> (SplitKKTMatrix::*)()>(&SplitKKTMatrix::hq))
    .def_property("hv", static_cast<const Eigen::VectorBlock<const Eigen::VectorXd> (SplitKKTMatrix::*)() const>(&SplitKKTMatrix::hv),
                        static_cast<Eigen::VectorBlock<Eigen::VectorXd> (SplitKKTMatrix::*)()>(&SplitKKTMatrix::hv))
    .def_readwrite("hu", &SplitKKTMatrix::hu)
    .def_readwrite("ha", &SplitKKTMatrix::ha)
    .def_property("hf", static_cast<const Eigen::VectorBlock<const Eigen::VectorXd> (SplitKKTMatrix::*)() const>(&SplitKKTMatrix::hf),
                        static_cast<Eigen::VectorBlock<Eigen::VectorXd> (SplitKKTMatrix::*)()>(&SplitKKTMatrix::hf))
    .def_property("Qff", static_cast<const Eigen::Block<const Eigen::MatrixXd> (SplitKKTMatrix::*)() const>(&SplitKKTMatrix::Qff),
                         static_cast<Eigen::Block<Eigen::MatrixXd> (SplitKKTMatrix::*)()>(&SplitKKTMatrix::Qff))
    .def_property("Qqf", static_cast<const Eigen::Block<const Eigen::MatrixXd> (SplitKKTMatrix::*)() const>(&SplitKKTMatrix::Qqf),
                         static_cast<Eigen::Block<Eigen::MatrixXd> (SplitKKTMatrix::*)()>(&SplitKKTMatrix::Qqf))
    .def("dimf", &SplitKKTMatrix::dimf)
    .def("dims", &SplitKKTMatrix::dims)
    DEFINE_ROBOTOC_PYBIND11_CLASS_CLONE(SplitKKTMatrix)
    DEFINE_ROBOTOC_PYBIND11_CLASS_PRINT(SplitKKTMatrix);
}

} // namespace python
} // namespace robotoc