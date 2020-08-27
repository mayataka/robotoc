#ifndef IDOCP_RICCATI_FACTORIZATION_HPP_
#define IDOCP_RICCATI_FACTORIZATION_HPP_

#include "Eigen/Core"

#include "idocp/robot/robot.hpp"


namespace idocp {

class RiccatiFactorization {
public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

  RiccatiFactorization(const Robot& robot) 
    : Pqq(Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv())),
      Pqv(Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv())),
      Pvq(Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv())),
      Pvv(Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv())),
      sq(Eigen::VectorXd::Zero(robot.dimv())),
      sv(Eigen::VectorXd::Zero(robot.dimv())) {
  }

  RiccatiFactorization() 
    : Pqq(),
      Pqv(),
      Pvq(),
      Pvv(),
      sq(),
      sv() {
  }

  ~RiccatiFactorization() {
  }

  RiccatiFactorization(const RiccatiFactorization&) = default;

  RiccatiFactorization& operator=(const RiccatiFactorization&) = default;
 
  RiccatiFactorization(RiccatiFactorization&&) noexcept = default;

  RiccatiFactorization& operator=(RiccatiFactorization&&) noexcept = default;

  Eigen::MatrixXd Pqq, Pqv, Pvq, Pvv;
  Eigen::VectorXd sq, sv;
};

} // namespace idocp 


#endif // IDOCP_RICCATI_FACTORIZATION_HPP_ 