#ifndef IDOCP_RICCATI_FACTORIZATION_HPP_
#define IDOCP_RICCATI_FACTORIZATION_HPP_

#include "Eigen/Core"

#include "idocp/robot/robot.hpp"


namespace idocp {

///
/// @class RiccatiFactorization
/// @brief Riccati factorized matrix and vector split into a time stage.
///
class RiccatiFactorization {
public:
  ///
  /// @brief Allocate Riccati factorization matrix and vector.
  /// @param[in] robot Robot model. Must be initialized by URDF or XML.
  ///
  RiccatiFactorization(const Robot& robot);

  ///
  /// @brief Default constructor. 
  ///
  RiccatiFactorization();

  ///
  /// @brief Destructor. 
  ///
  ~RiccatiFactorization();

  ///
  /// @brief Default copy constructor. 
  ///
  RiccatiFactorization(const RiccatiFactorization&) = default;

  ///
  /// @brief Default copy operator. 
  ///
  RiccatiFactorization& operator=(const RiccatiFactorization&) = default;

  ///
  /// @brief Default move constructor. 
  ///
  RiccatiFactorization(RiccatiFactorization&&) noexcept = default;

  ///
  /// @brief Default move assign operator. 
  ///
  RiccatiFactorization& operator=(RiccatiFactorization&&) noexcept = default;

  ///
  /// @brief Riccati factorization matrix. Size is 
  /// Robot::dimv() x Robot::dimv().
  ///
  Eigen::MatrixXd Pqq;

  ///
  /// @brief Riccati factorization matrix. Size is 
  /// Robot::dimv() x Robot::dimv().
  ///
  Eigen::MatrixXd Pqv;

  ///
  /// @brief Riccati factorization matrix. Size is 
  /// Robot::dimv() x Robot::dimv().
  ///
  Eigen::MatrixXd Pvq;

  ///
  /// @brief Riccati factorization matrix. Size is 
  /// Robot::dimv() x Robot::dimv().
  ///
  Eigen::MatrixXd Pvv;

  ///
  /// @brief Riccati factorization vector. Size is Robot::dimv().
  ///
  Eigen::VectorXd sq;

  ///
  /// @brief Riccati factorization vector. Size is Robot::dimv().
  ///
  Eigen::VectorXd sv;

  ///
  /// @brief Feedback gain of LQR subproblem. Size is 
  /// Robot::dimu() x 2 * Robot::dimv().
  ///
  Eigen::MatrixXd K;

  ///
  /// @brief Feedforward term of LQR subproblem. Size is Robot::dimu().
  ///
  Eigen::VectorXd k;

  ///
  /// @brief Riccati factorization for pure-state equality constraints. 
  /// Size is 2 * Robot::dimv() x 2 * Robot::dimv() if the robot can have
  /// contacts. Size is 0 x 0 otherwise.
  ///
  Eigen::MatrixXd ApBK;

  ///
  /// @brief Riccati factorization for pure-state equality constraints. Size is 
  /// 2 * Robot::dimv(). Size is 0 otherwise.
  ///
  Eigen::VectorXd apBk;

  ///
  /// @brief Riccati factorization for pure-state equality constraints. 
  /// Size is Robot::dimv() x Robot::dimv() if the robot can have
  /// contacts. Size is 0 x 0 otherwise.
  ///
  Eigen::MatrixXd BGinvBt;

  ///
  /// @brief Riccati factorization for pure-state equality constraints. 
  /// Size is 2 * Robot::dimv() x 2 * Robot::dimv() if the robot can have
  /// contacts. Size is 0 x 0 otherwise.
  ///
  Eigen::MatrixXd Pi;

  ///
  /// @brief Riccati factorization for pure-state equality constraints. 
  /// Size is 2 * Robot::dimv() x 2 * Robot::dimv() if the robot can have
  /// contacts. Size is 0 x 0 otherwise.
  ///
  Eigen::MatrixXd N;

  ///
  /// @brief Riccati factorization for pure-state equality constraints. Size is 
  /// 2 * Robot::dimv(). Size is 0 otherwise.
  ///
  Eigen::VectorXd pi;

  ///
  /// @brief Left cols of the feedback gain of the LQR subproblem. Size is 
  /// Robot::dimu() x Robot::dimv().
  ///
  Eigen::Block<Eigen::MatrixXd> Kq();

  ///
  /// @brief const version of RiccatiFactorization::Kq() 
  ///
  const Eigen::Block<const Eigen::MatrixXd> Kq() const;

  ///
  /// @brief Right cols of the feedback gain of the LQR subproblem. Size is 
  /// Robot::dimu() x Robot::dimv().
  ///
  Eigen::Block<Eigen::MatrixXd> Kv();

  ///
  /// @brief const version of RiccatiFactorization::Kv() 
  ///
  const Eigen::Block<const Eigen::MatrixXd> Kv() const;

  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

private:
  int dimv_, dimu_;

};

} // namespace idocp 

#include "idocp/ocp/riccati_factorization.hxx"

#endif // IDOCP_RICCATI_FACTORIZATION_HPP_ 