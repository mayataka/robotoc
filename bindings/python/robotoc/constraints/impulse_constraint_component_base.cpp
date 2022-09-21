#include <pybind11/pybind11.h>

#include "robotoc/constraints/impulse_constraint_component_base.hpp"


namespace robotoc {
namespace python {

class PyImpulseConstraintComponentBase : public ImpulseConstraintComponentBase {
public:
  // Inherit the constructors
  using ImpulseConstraintComponentBase::ImpulseConstraintComponentBase;

  KinematicsLevel kinematicsLevel() const override {
    PYBIND11_OVERRIDE_PURE(KinematicsLevel, ImpulseConstraintComponentBase, 
                           kinematicsLevel, );
  }

  void allocateExtraData(ConstraintComponentData& data) const override {
    PYBIND11_OVERRIDE_PURE(void, ImpulseConstraintComponentBase, 
                           allocateExtraData, data);
  }

  bool isFeasible(Robot& robot, const ImpulseStatus& impulse_status, 
                  ConstraintComponentData& data, 
                  const SplitSolution& s) const override {
    PYBIND11_OVERRIDE_PURE(bool, ImpulseConstraintComponentBase, 
                           isFeasible, 
                           robot, impulse_status, data, s);
  }

  void setSlack(Robot& robot, const ImpulseStatus& impulse_status, 
                ConstraintComponentData& data, 
                const SplitSolution& s) const override {
    PYBIND11_OVERRIDE_PURE(void, ImpulseConstraintComponentBase, 
                           setSlack, 
                           robot, impulse_status, data, s);
  }

  void evalConstraint(Robot& robot, const ImpulseStatus& impulse_status, 
                      ConstraintComponentData& data, 
                      const SplitSolution& s) const override {
    PYBIND11_OVERRIDE_PURE(void, ImpulseConstraintComponentBase, 
                           evalConstraint, 
                           robot, impulse_status, data, s);
  }

  void evalDerivatives(Robot& robot, const ImpulseStatus& impulse_status, 
                       ConstraintComponentData& data,
                       const SplitSolution& s,
                       SplitKKTResidual& kkt_residual) const override {
    PYBIND11_OVERRIDE_PURE(void, ImpulseConstraintComponentBase, 
                           evalDerivatives, 
                           robot, impulse_status, data, s, kkt_residual);
  }

  void condenseSlackAndDual(const ImpulseStatus& impulse_status, 
                            ConstraintComponentData& data,
                            SplitKKTMatrix& kkt_matrix,
                            SplitKKTResidual& kkt_residual) const override {
    PYBIND11_OVERRIDE_PURE(void, ImpulseConstraintComponentBase, 
                           condenseSlackAndDual, 
                           impulse_status, data, kkt_matrix, kkt_residual);
  }

  void expandSlackAndDual(const ImpulseStatus& impulse_status, 
                          ConstraintComponentData& data, 
                          const SplitDirection& d) const override {
    PYBIND11_OVERRIDE_PURE(void, ImpulseConstraintComponentBase, 
                           expandSlackAndDual, 
                           impulse_status, data, d);
  }

  int dimc() const override {
    PYBIND11_OVERRIDE_PURE(int, ImpulseConstraintComponentBase, 
                           dimc, );
  }
};


namespace py = pybind11;

PYBIND11_MODULE(impulse_constraint_component_base, m) {
  py::class_<ImpulseConstraintComponentBase, 
             PyImpulseConstraintComponentBase, 
             std::shared_ptr<ImpulseConstraintComponentBase>>(m, "ImpulseConstraintComponentBase")
    .def(py::init<const double, const double>(),
          py::arg("barrier_param")=1.0e-03, py::arg("fraction_to_boundary_rule")=0.995)
    .def("kinematicsLevel", &ImpulseConstraintComponentBase::kinematicsLevel)
    .def("allocateExtraData", &ImpulseConstraintComponentBase::allocateExtraData,
          py::arg("data"))
    .def("isFeasible", &ImpulseConstraintComponentBase::isFeasible,
          py::arg("robot"), py::arg("contact_status"), py::arg("data"), py::arg("s"))
    .def("setSlack", &ImpulseConstraintComponentBase::setSlack,
          py::arg("robot"), py::arg("contact_status"), py::arg("data"), py::arg("s"))
    .def("evalConstraint", &ImpulseConstraintComponentBase::evalConstraint,
          py::arg("robot"), py::arg("contact_status"), py::arg("data"), py::arg("s"))
    .def("evalDerivatives", &ImpulseConstraintComponentBase::evalDerivatives,
          py::arg("robot"), py::arg("contact_status"), py::arg("data"), py::arg("s"),
          py::arg("kkt_residual"))
    .def("condenseSlackAndDual", &ImpulseConstraintComponentBase::condenseSlackAndDual,
          py::arg("contact_status"), py::arg("data"), py::arg("kkt_matrix"), 
          py::arg("kkt_residual"))
    .def("expandSlackAndDual", &ImpulseConstraintComponentBase::expandSlackAndDual,
          py::arg("contact_status"), py::arg("data"), py::arg("d"))
    .def("dimc", &ImpulseConstraintComponentBase::dimc);
}

} // namespace python
} // namespace robotoc