#ifndef ROBOTOC_STATE_EQUATION_DATA_HPP_
#define ROBOTOC_STATE_EQUATION_DATA_HPP_

#include "Eigen/Core"

#include "robotoc/robot/robot.hpp"
#include "robotoc/robot/se3_jacobian_inverse.hpp"


namespace robotoc {

///
/// @class StateEquationData
/// @brief Data for the state equations. 
///
class StateEquationData {
public:
  ///
  /// @brief Constructs a state equation.
  /// @param[in] robot Robot model. 
  ///
  StateEquationData(const Robot& robot);

  ///
  /// @brief Default constructor.  
  ///
  StateEquationData();
  
  ///
  /// @brief Destructor. 
  ///
  ~StateEquationData() = default;

  ///
  /// @brief Default copy constructor. 
  ///
  StateEquationData(const StateEquationData&) = default;

  ///
  /// @brief Default copy assign operator. 
  ///
  StateEquationData& operator=(const StateEquationData&) = default;

  ///
  /// @brief Default move constructor. 
  ///
  StateEquationData(StateEquationData&&) noexcept = default;

  ///
  /// @brief Default move assign operator. 
  ///
  StateEquationData& operator=(StateEquationData&&) noexcept = default;

  Eigen::MatrixXd Fqq_prev;

  Eigen::MatrixXd Fqq_inv;

  Eigen::MatrixXd Fqq_prev_inv;

  Eigen::MatrixXd Fqq_tmp;  

  Eigen::VectorXd Fq_tmp;

  SE3JacobianInverse se3_jac_inverse;

  bool has_floating_base;
};

} // namespace robotoc 

#endif // ROBOTOC_STATE_EQUATION_DATA_HPP_