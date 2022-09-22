#include "robotoc/dynamics/state_equation_data.hpp"


namespace robotoc {

StateEquationData::StateEquationData(const Robot& robot) 
  : Fqq_prev(Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv())),
    Fqq_inv(),
    Fqq_prev_inv(),
    Fqq_tmp(),
    Fq_tmp(),
    se3_jac_inverse(),
    has_floating_base(robot.hasFloatingBase()) {
  if (robot.hasFloatingBase()) {
    Fqq_inv_.resize(6, 6);
    Fqq_inv_.setZero();
    Fqq_prev_inv_.resize(6, 6);
    Fqq_prev_inv_.setZero();
    Fqq_tmp_.resize(6, 6);
    Fqq_tmp_.setZero();
    Fq_tmp_.resize(6);
    Fq_tmp_.setZero();
  }
}


StateEquationData::StateEquationData() 
  : Fqq_prev(),
    Fqq_inv(),
    Fqq_prev_inv(),
    Fqq_tmp(),
    Fq_tmp(),
    se3_jac_inverse(),
    has_floating_base(false) {
}

} // namespace robotoc 