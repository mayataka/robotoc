#include <pybind11/pybind11.h>

#include "idocp/cost/cost_function_component_base.hpp"


namespace idocp {
namespace python {

class PyCostFunctionComponentBase : public CostFunctionComponentBase {
public:
  // Inherit the constructors
  using CostFunctionComponentBase::CostFunctionComponentBase;

  bool useKinematics() const override {
    PYBIND11_OVERRIDE_PURE(bool, CostFunctionComponentBase, 
                           useKinematics, );
  }

  double computeStageCost(Robot& robot, CostFunctionData& data, const double t, 
                          const double dt, 
                          const SplitSolution& s) const override {
    PYBIND11_OVERRIDE_PURE(double, CostFunctionComponentBase, 
                           computeStageCost, 
                           robot, data, t, dt, s);
  }

  void computeStageCostDerivatives(
      Robot& robot, CostFunctionData& data, const double t, const double dt, 
      const SplitSolution& s, SplitKKTResidual& kkt_residual) const override {
    PYBIND11_OVERRIDE_PURE(void, CostFunctionComponentBase, 
                           computeStageCostDerivatives, 
                           robot, data, t, dt, s, kkt_residual);
  }

  void computeStageCostHessian(
      Robot& robot, CostFunctionData& data, const double t, const double dt, 
      const SplitSolution& s, SplitKKTMatrix& kkt_matrix) const override {
    PYBIND11_OVERRIDE_PURE(void, CostFunctionComponentBase, 
                           computeStageCostHessian, 
                           robot, data, t, dt, s, kkt_matrix);
  }

  double computeTerminalCost(Robot& robot, CostFunctionData& data, 
                             const double t, 
                             const SplitSolution& s) const override {
    PYBIND11_OVERRIDE_PURE(double, CostFunctionComponentBase, 
                           computeTerminalCost, 
                           robot, data, t, s);
  }

  void computeTerminalCostDerivatives(
      Robot& robot, CostFunctionData& data, const double t, 
      const SplitSolution& s, SplitKKTResidual& kkt_residual) const override {
    PYBIND11_OVERRIDE_PURE(void, CostFunctionComponentBase, 
                           computeTerminalCostDerivatives, 
                           robot, data, t, s, kkt_residual);
  }

  void computeTerminalCostHessian(
      Robot& robot, CostFunctionData& data, const double t, 
      const SplitSolution& s, SplitKKTMatrix& kkt_matrix) const override {
    PYBIND11_OVERRIDE_PURE(void, CostFunctionComponentBase, 
                           computeTerminalCostHessian, 
                           robot, data, t, s, kkt_matrix);
  }

  double computeImpulseCost(Robot& robot, CostFunctionData& data, 
                            const double t, 
                            const ImpulseSplitSolution& s) const override {
    PYBIND11_OVERRIDE_PURE(double, CostFunctionComponentBase, 
                           computeImpulseCost, 
                           robot, data, t, s);
  }

  void computeImpulseCostDerivatives(
      Robot& robot, CostFunctionData& data, const double t, 
      const ImpulseSplitSolution& s, 
      ImpulseSplitKKTResidual& kkt_residual) const override {
    PYBIND11_OVERRIDE_PURE(void, CostFunctionComponentBase, 
                           computeImpulseCostDerivatives, 
                           robot, data, t, s, kkt_residual);
  }

  void computeImpulseCostHessian(
      Robot& robot, CostFunctionData& data, const double t, 
      const ImpulseSplitSolution& s, 
      ImpulseSplitKKTMatrix& kkt_matrix) const override {
    PYBIND11_OVERRIDE_PURE(void, CostFunctionComponentBase, 
                           computeImpulseCostHessian, 
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
} // namespace idocp