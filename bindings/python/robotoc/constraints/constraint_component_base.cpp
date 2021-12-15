#include <pybind11/pybind11.h>

#include "robotoc/constraints/constraint_component_base.hpp"


namespace robotoc {
namespace python {

class PyConstraintComponentBase : public ConstraintComponentBase {
public:
  // Inherit the constructors
  using ConstraintComponentBase::ConstraintComponentBase;

  bool useKinematics() const override {
    PYBIND11_OVERRIDE_PURE(bool, ConstraintComponentBase, useKinematics, );
  }

  KinematicsLevel kinematicsLevel() const override {
    PYBIND11_OVERRIDE_PURE(KinematicsLevel, ConstraintComponentBase, 
                           kinematicsLevel, );
  }

  void allocateExtraData(ConstraintComponentData& data) const override {
    PYBIND11_OVERRIDE_PURE(void, ConstraintComponentBase, 
                           allocateExtraData, 
                           data);
  }

  bool isFeasible(Robot& robot, const ContactStatus& contact_status, 
                  ConstraintComponentData& data, 
                  const SplitSolution& s) const override {
    PYBIND11_OVERRIDE_PURE(bool, ConstraintComponentBase, 
                           isFeasible, 
                           robot, contact_status, data, s);
  }

  void setSlack(Robot& robot, const ContactStatus& contact_status, 
                ConstraintComponentData& data, 
                const SplitSolution& s) const override {
    PYBIND11_OVERRIDE_PURE(void, ConstraintComponentBase, 
                           setSlack, 
                           robot, contact_status, data, s);
  }

  void evalConstraint(Robot& robot, const ContactStatus& contact_status, 
                      ConstraintComponentData& data, 
                      const SplitSolution& s) const override {
    PYBIND11_OVERRIDE_PURE(void, ConstraintComponentBase, 
                           evalConstraint, 
                           robot, contact_status, data, s);
  }

  void evalDerivatives(Robot& robot, const ContactStatus& contact_status, 
                       ConstraintComponentData& data, const SplitSolution& s,
                       SplitKKTResidual& kkt_residual) const override {
    PYBIND11_OVERRIDE_PURE(void, ConstraintComponentBase, 
                           evalDerivatives, 
                           robot, contact_status, data, s, kkt_residual);
  }

  void condenseSlackAndDual(const ContactStatus& contact_status,
                            ConstraintComponentData& data, 
                            SplitKKTMatrix& kkt_matrix,
                            SplitKKTResidual& kkt_residual) const override {
    PYBIND11_OVERRIDE_PURE(void, ConstraintComponentBase, 
                           condenseSlackAndDual, 
                           contact_status, data, kkt_matrix, kkt_residual);
  }

  void expandSlackAndDual(const ContactStatus& contact_status, 
                          ConstraintComponentData& data, 
                          const SplitDirection& d) const override {
    PYBIND11_OVERRIDE_PURE(void, ConstraintComponentBase, 
                           expandSlackAndDual, 
                           contact_status, data, d);
  }

  int dimc() const override {
    PYBIND11_OVERRIDE_PURE(int, ConstraintComponentBase, 
                           dimc, );
  }
};


namespace py = pybind11;

PYBIND11_MODULE(constraint_component_base, m) {
  py::enum_<KinematicsLevel>(m, "KinematicsLevel", py::arithmetic())
    .value("PositionLevel",  KinematicsLevel::PositionLevel)
    .value("VelocityLevel", KinematicsLevel::VelocityLevel)
    .value("AccelerationLevel", KinematicsLevel::AccelerationLevel)
    .export_values();

  py::class_<ConstraintComponentBase, 
             PyConstraintComponentBase,
             std::shared_ptr<ConstraintComponentBase>>(m, "ConstraintComponentBase")
    .def(py::init<const double, const double>(),
          py::arg("barrier")=1.0e-03, py::arg("fraction_to_boundary_rule")=0.995);
}

} // namespace python
} // namespace robotoc