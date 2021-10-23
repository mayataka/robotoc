#include <pybind11/pybind11.h>

#include "robotoc/hybrid/switching_time_cost_function_component_base.hpp"


namespace robotoc {
namespace python {

class PySwitchingTimeCostFunctionComponentBase : public SwitchingTimeCostFunctionComponentBase {
public:
  // Inherit the constructors
  using SwitchingTimeCostFunctionComponentBase::SwitchingTimeCostFunctionComponentBase;

  double computeCost(const double t0, const double tf, 
                     const Eigen::VectorXd& ts) const override {
    PYBIND11_OVERRIDE_PURE(double, SwitchingTimeCostFunctionComponentBase, 
                           computeCost, 
                           t0, tf, ts);
  }

  void computeCostDerivatives(const double t0, const double tf, 
                              const Eigen::VectorXd& ts,
                              Eigen::VectorXd& hts) const override {
    PYBIND11_OVERRIDE_PURE(void, SwitchingTimeCostFunctionComponentBase, 
                           computeCostDerivatives, 
                           t0, tf, ts, hts);
  }

  void computeCostHessian(const double t0, const double tf, 
                          const Eigen::VectorXd& ts,
                          Eigen::MatrixXd& Qts) const override {
    PYBIND11_OVERRIDE_PURE(void, SwitchingTimeCostFunctionComponentBase, 
                           computeStageCostHessian, 
                           t0, tf, ts, Qts);
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