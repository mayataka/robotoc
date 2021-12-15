#include <pybind11/pybind11.h>

#include "robotoc/cost/cost_function_component_base.hpp"


namespace robotoc {
namespace python {

class PyCostFunctionComponentBase : public CostFunctionComponentBase {
public:
  // Inherit the constructors
  using CostFunctionComponentBase::CostFunctionComponentBase;

  bool useKinematics() const override {
    PYBIND11_OVERRIDE_PURE(bool, CostFunctionComponentBase, 
                           useKinematics, );
  }

  double evalStageCost(Robot& robot, const ContactStatus& contact_status, 
                       CostFunctionData& data, const double t, const double dt, 
                       const SplitSolution& s) const override {
    PYBIND11_OVERRIDE_PURE(double, CostFunctionComponentBase, 
                           evalStageCost, 
                           robot, contact_status, data, t, dt, s);
  }

  void evalStageCostDerivatives(Robot& robot, const ContactStatus& contact_status, 
                                CostFunctionData& data, const double t, 
                                const double dt, const SplitSolution& s, 
                                SplitKKTResidual& kkt_residual) const override {
    PYBIND11_OVERRIDE_PURE(void, CostFunctionComponentBase, 
                           evalStageCostDerivatives, 
                           robot, contact_status, data, t, dt, s, kkt_residual);
  }

  void evalStageCostHessian(Robot& robot, const ContactStatus& contact_status, 
                            CostFunctionData& data, const double t, 
                            const double dt, const SplitSolution& s, 
                            SplitKKTMatrix& kkt_matrix) const override {
    PYBIND11_OVERRIDE_PURE(void, CostFunctionComponentBase, 
                           evalStageCostHessian, 
                           robot, contact_status, data, t, dt, s, kkt_matrix);
  }

  double evalTerminalCost(Robot& robot, CostFunctionData& data, 
                          const double t, 
                          const SplitSolution& s) const override {
    PYBIND11_OVERRIDE_PURE(double, CostFunctionComponentBase, 
                           evalTerminalCost, 
                           robot, data, t, s);
  }

  void evalTerminalCostDerivatives(Robot& robot, CostFunctionData& data, 
                                   const double t, const SplitSolution& s, 
                                   SplitKKTResidual& kkt_residual) const override {
    PYBIND11_OVERRIDE_PURE(void, CostFunctionComponentBase, 
                           evalTerminalCostDerivatives, 
                           robot, data, t, s, kkt_residual);
  }

  void evalTerminalCostHessian(Robot& robot, CostFunctionData& data, 
                               const double t, const SplitSolution& s, 
                               SplitKKTMatrix& kkt_matrix) const override {
    PYBIND11_OVERRIDE_PURE(void, CostFunctionComponentBase, 
                           evalTerminalCostHessian, 
                           robot, data, t, s, kkt_matrix);
  }

  double evalImpulseCost(Robot& robot, const ImpulseStatus& impulse_status, 
                         CostFunctionData& data, const double t, 
                         const ImpulseSplitSolution& s) const override {
    PYBIND11_OVERRIDE_PURE(double, CostFunctionComponentBase, 
                           evalImpulseCost, 
                           robot, impulse_status, data, t, s);
  }

  void evalImpulseCostDerivatives(Robot& robot, const ImpulseStatus& impulse_status, 
                                  CostFunctionData& data, const double t, 
                                  const ImpulseSplitSolution& s, 
                                  ImpulseSplitKKTResidual& kkt_residual) const override {
    PYBIND11_OVERRIDE_PURE(void, CostFunctionComponentBase, 
                           evalImpulseCostDerivatives, 
                           robot, impulse_status, data, t, s, kkt_residual);
  }

  void evalImpulseCostHessian(Robot& robot, const ImpulseStatus& impulse_status, 
                              CostFunctionData& data, const double t, 
                              const ImpulseSplitSolution& s, 
                              ImpulseSplitKKTMatrix& kkt_matrix) const override {
    PYBIND11_OVERRIDE_PURE(void, CostFunctionComponentBase, 
                           evalImpulseCostHessian, 
                           robot, impulse_status, data, t, s, kkt_matrix);
  }
};


namespace py = pybind11;

PYBIND11_MODULE(cost_function_component_base, m) {
  py::class_<CostFunctionComponentBase, 
             PyCostFunctionComponentBase, 
             std::shared_ptr<CostFunctionComponentBase>>(m, "CostFunctionComponentBase")
    .def(py::init<>());

}

} // namespace python
} // namespace robotoc