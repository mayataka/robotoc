#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include <pybind11/numpy.h>

#include "robotoc/hybrid/sto_cost_function_component_base.hpp"
#include "robotoc/hybrid/time_discretization.hpp"


namespace robotoc {
namespace python {

class PySTOCostFunctionComponentBase : public STOCostFunctionComponentBase {
public:
  // Inherit the constructors
  using STOCostFunctionComponentBase::STOCostFunctionComponentBase;

  double evalCost(const TimeDiscretization& discretization) const override {
    PYBIND11_OVERRIDE_PURE(double, STOCostFunctionComponentBase, 
                           evalCost, 
                           discretization);
  }

  void evalCostDerivatives(const TimeDiscretization& discretization,
                           Eigen::VectorXd& lts) const override {
    PYBIND11_OVERRIDE_PURE(void, STOCostFunctionComponentBase, 
                           evalCostDerivatives, 
                           discretization, lts);
  }

  void evalCostHessian(const TimeDiscretization& discretization,
                       Eigen::MatrixXd& Qts) const override {
    PYBIND11_OVERRIDE_PURE(void, STOCostFunctionComponentBase, 
                           evalStageCostHessian, 
                           discretization, Qts);
  }
};


namespace py = pybind11;

PYBIND11_MODULE(sto_cost_function_component_base, m) {
  py::class_<STOCostFunctionComponentBase, 
             PySTOCostFunctionComponentBase, 
             std::shared_ptr<STOCostFunctionComponentBase>>(m, "STOCostFunctionComponentBase")
    .def(py::init<>());

}

} // namespace python
} // namespace robotoc