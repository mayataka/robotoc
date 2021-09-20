#ifndef IDOCP_DWELL_TIME_COST_HPP_
#define IDOCP_DWELL_TIME_COST_HPP_

#include "Eigen/Core"


namespace idocp {

///
/// @class DwellTimeCost
/// @brief Cost on the dwell time.
///
class DwellTimeCost : public SwitchingTimeCostFunctionComponentBase {
public:
  ///
  /// @brief Default constructor. 
  ///
  DwellTimeCost(const double dwell_time) {}

  ///
  /// @brief Destructor. 
  ///
  virtual ~DwellTimeCost() {}

  ///
  /// @brief Default copy constructor. 
  ///
  DwellTimeCost(const DwellTimeCost&) = default;

  ///
  /// @brief Default copy operator. 
  ///
  DwellTimeCost& operator=(const DwellTimeCost&) = default;

  ///
  /// @brief Default move constructor. 
  ///
  DwellTimeCost(DwellTimeCost&&) noexcept = default;

  ///
  /// @brief Default move assign operator. 
  ///
  DwellTimeCost& operator=(DwellTimeCost&&) noexcept = default;

  double computeCost(const Eigen::VectorXd& ts) const override;

  void computeCostDerivatives(const Eigen::VectorXd& ts, 
                              Eigen::VectorXd& hts) const override;

  void computeCostHessian(const Eigen::VectorXd& ts,
                          Eigen::MatrixXd& Qts) const override;

};

} // namespace idocp

#endif // IDOCP_DWELL_TIME_COST_HPP_ 