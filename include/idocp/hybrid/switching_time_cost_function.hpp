#ifndef IDOCP_SWITCHING_TIME_COST_FUNCTION_HPP_
#define IDOCP_SWITCHING_TIME_COST_FUNCTION_HPP_

#include <vector>
#include <memory>

#include "Eigen/Core"

#include "idocp/hybrid/hybrid_time_discretization.hpp"
#include "idocp/hybrid/switching_time_cost_function_component_base.hpp"


namespace idocp {

///
/// @class SwitchingTimeCostFunction
/// @brief Stack of the cost function on the switchig times. Composed by cost 
/// function components that inherits SwitchingTimeCostFunctionComponentBase.
///
class SwitchingTimeCostFunction {
public:
  using SwitchingTimeCostFunctionComponentBasePtr 
      = std::shared_ptr<SwitchingTimeCostFunctionComponentBase>;

  ///
  /// @brief Default constructor. 
  ///
  SwitchingTimeCostFunction();

  ///
  /// @brief Destructor. 
  ///
  ~SwitchingTimeCostFunction();

  ///
  /// @brief Default copy constructor. 
  ///
  SwitchingTimeCostFunction(const SwitchingTimeCostFunction&) = default;

  ///
  /// @brief Default copy operator. 
  ///
  SwitchingTimeCostFunction& operator=(
      const SwitchingTimeCostFunction&) = default;

  ///
  /// @brief Default move constructor. 
  ///
  SwitchingTimeCostFunction(SwitchingTimeCostFunction&&) noexcept = default;

  ///
  /// @brief Default move assign operator. 
  ///
  SwitchingTimeCostFunction& operator=(
      SwitchingTimeCostFunction&&) noexcept = default;

  ///
  /// @brief Append a cost function component to the cost function.
  /// @param[in] cost shared pointer to the switching tmie cost function  
  /// component appended to the cost.
  ///
  void push_back(const SwitchingTimeCostFunctionComponentBasePtr& cost);

  ///
  /// @brief Clear cost function by removing all components.
  ///
  void clear();

  ///
  /// @brief Computes the cost on the switching times. 
  /// @param[in] discretization Discretization of the optimal control problem.
  /// @return Cost on the switching times.
  ///
  double computeCost(const HybridTimeDiscretization& discretization);

  ///
  /// @brief Computes the cost on the switching times and its first-order 
  /// partial derivatives. 
  /// @param[in] discretization Discretization of the optimal control problem.
  /// @param[in, out] kkt_residual KKT residual. The partial derivatives 
  /// are added to this object.
  /// @return Cost on the switching times.
  ///
  double linearizeCost(const HybridTimeDiscretization& discretization,
                       KKTResidual& kkt_residual); 

  ///
  /// @brief Computes the cost, its first-order partial derivatives, and 
  /// its Hessian, i.e., its second-order partial derivatives. 
  /// @param[in] discretization Discretization of the optimal control problem.
  /// @param[out] kkt_matrix KKT matrix.
  /// @param[out] kkt_residual KKT residual.
  /// @return Cost on the switching times.
  ///
  double quadratizeStageCost(const HybridTimeDiscretization& discretization,
                             KKTMatrix& kkt_matrix, 
                             KKTResidual& kkt_residual);

private:
  std::vector<SwitchingTimeCostFunctionComponentBasePtr> costs_;
  Eigen::VectorXd ts_, hts_;
  Eigen::MatrixXd Qts_;

  void setNumSwitches(const int num_switches);

  void setSwitchingTimes(const HybridTimeDiscretization& discretization);

  void setKKT(const HybridTimeDiscretization& discretization,
              KKTMatrix& kkt_matrix, KKTResidual& kkt_residual);

  void setKKT(const HybridTimeDiscretization& discretization,
              KKTResidual& kkt_residual);

};

} // namespace idocp

#include "idocp/hybrid/switching_time_cost_function.hxx"

#endif // IIDOCP_SWITCHING_TIME_COST_FUNCTION_HPP