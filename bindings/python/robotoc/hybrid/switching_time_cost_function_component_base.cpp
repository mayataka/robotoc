#include <pybind11/pybind11.h>

#include "robotoc/hybrid/switching_time_cost_function_component_base.hpp"
#include "robotoc/hybrid/hybrid_ocp_discretization.hpp"


namespace robotoc {
namespace python {

class PySwitchingTimeCostFunctionComponentBase : public SwitchingTimeCostFunctionComponentBase {
public:
  // Inherit the constructors
  using SwitchingTimeCostFunctionComponentBase::SwitchingTimeCostFunctionComponentBase;

  double evalCost(const HybridOCPDiscretization& discretization) const override {
    PYBIND11_OVERRIDE_PURE(double, SwitchingTimeCostFunctionComponentBase, 
                           evalCost, 
                           discretization);
  }

  void evalCostDerivatives(const HybridOCPDiscretization& discretization,
                           Eigen::VectorXd& lts) const override {
    PYBIND11_OVERRIDE_PURE(void, SwitchingTimeCostFunctionComponentBase, 
                           evalCostDerivatives, 
                           discretization, lts);
  }

  void evalCostHessian(const HybridOCPDiscretization& discretization,
                       Eigen::MatrixXd& Qts) const override {
    PYBIND11_OVERRIDE_PURE(void, SwitchingTimeCostFunctionComponentBase, 
                           evalStageCostHessian, 
                           discretization, Qts);
  }
};


namespace py = pybind11;

PYBIND11_MODULE(switching_time_cost_function_component_base, m) {
  py::class_<SwitchingTimeCostFunctionComponentBase, 
             PySwitchingTimeCostFunctionComponentBase, 
             std::shared_ptr<SwitchingTimeCostFunctionComponentBase>>(m, "SwitchingTimeCostFunctionComponentBase")
    .def(py::init<>());

}

} // namespace python
} // namespace robotoc