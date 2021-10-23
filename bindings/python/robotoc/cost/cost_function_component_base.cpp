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

  double evalStageCost(Robot& robot, CostFunctionData& data, const double t, 
                       const double dt, const SplitSolution& s) const override {
    PYBIND11_OVERRIDE_PURE(double, CostFunctionComponentBase, 
                           evalStageCost, 
                           robot, data, t, dt, s);
  }

  void evalStageCostDerivatives(Robot& robot, CostFunctionData& data, 
                                const double t, const double dt, 
                                const SplitSolution& s, 
                                SplitKKTResidual& kkt_residual) const override {
    PYBIND11_OVERRIDE_PURE(void, CostFunctionComponentBase, 
                           evalStageCostDerivatives, 
                           robot, data, t, dt, s, kkt_residual);
  }

  void evalStageCostHessian(Robot& robot, CostFunctionData& data, 
                            const double t, const double dt, 
                            const SplitSolution& s, 
                            SplitKKTMatrix& kkt_matrix) const override {
    PYBIND11_OVERRIDE_PURE(void, CostFunctionComponentBase, 
                           evalStageCostHessian, 
                           robot, data, t, dt, s, kkt_matrix);
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

  double evalImpulseCost(Robot& robot, CostFunctionData& data, const double t, 
                         const ImpulseSplitSolution& s) const override {
    PYBIND11_OVERRIDE_PURE(double, CostFunctionComponentBase, 
                           evalImpulseCost, 
                           robot, data, t, s);
  }

  void evalImpulseCostDerivatives(Robot& robot, CostFunctionData& data, 
                                  const double t, const ImpulseSplitSolution& s, 
                                  ImpulseSplitKKTResidual& kkt_residual) const override {
    PYBIND11_OVERRIDE_PURE(void, CostFunctionComponentBase, 
                           evalImpulseCostDerivatives, 
                           robot, data, t, s, kkt_residual);
  }

  void evalImpulseCostHessian(Robot& robot, CostFunctionData& data, const double t, 
                              const ImpulseSplitSolution& s, 
                              ImpulseSplitKKTMatrix& kkt_matrix) const override {
    PYBIND11_OVERRIDE_PURE(void, CostFunctionComponentBase, 
                           evalImpulseCostHessian, 
                           robot, data, t, s, kkt_matrix);
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