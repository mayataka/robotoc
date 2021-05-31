#include <pybind11/pybind11.h>

#include "idocp/constraints/impulse_constraint_component_base.hpp"


namespace idocp {
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

  bool isFeasible(Robot& robot, ConstraintComponentData& data, 
                  const ImpulseSplitSolution& s) const override {
    PYBIND11_OVERRIDE_PURE(bool, ImpulseConstraintComponentBase, 
                           isFeasible, 
                           robot, data, s);
  }

  void setSlackAndDual(Robot& robot, ConstraintComponentData& data, 
                       const ImpulseSplitSolution& s) const override {
    PYBIND11_OVERRIDE_PURE(void, ImpulseConstraintComponentBase, 
                           setSlackAndDual, 
                           robot, data, s);
  }

  void augmentDualResidual(Robot& robot, ConstraintComponentData& data,
                           const ImpulseSplitSolution& s,
                           ImpulseSplitKKTResidual& kkt_residual) const override {
    PYBIND11_OVERRIDE_PURE(void, ImpulseConstraintComponentBase, 
                           augmentDualResidual, 
                           robot, data, s, kkt_residual);
  }

  void condenseSlackAndDual(Robot& robot, ConstraintComponentData& data,
                            const ImpulseSplitSolution& s, 
                            ImpulseSplitKKTMatrix& kkt_matrix,
                            ImpulseSplitKKTResidual& kkt_residual) const override {
    PYBIND11_OVERRIDE_PURE(void, ImpulseConstraintComponentBase, 
                           condenseSlackAndDual, 
                           robot, data, s, kkt_matrix, kkt_residual);
  }

  void expandSlackAndDual(ConstraintComponentData& data, 
                          const ImpulseSplitSolution& s, 
                          const ImpulseSplitDirection& d) const override {
    PYBIND11_OVERRIDE_PURE(void, ImpulseConstraintComponentBase, 
                           expandSlackAndDual, 
                           data, s, d);
  }

  void computePrimalAndDualResidual(Robot& robot, ConstraintComponentData& data, 
                                    const ImpulseSplitSolution& s) const override {
    PYBIND11_OVERRIDE_PURE(void, ImpulseConstraintComponentBase, 
                           computePrimalAndDualResidual, 
                           robot, data, s);
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
    .def(py::init<const double, const double>());
}

} // namespace python
} // namespace idocp