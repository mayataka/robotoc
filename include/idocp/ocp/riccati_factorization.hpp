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
  RiccatiFactorization(const Robot& robot)
    : Pqq(Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv())),
      Pqv(Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv())),
      Pvq(Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv())),
      Pvv(Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv())),
      sq(Eigen::VectorXd::Zero(robot.dimv())),
      sv(Eigen::VectorXd::Zero(robot.dimv())),
      Pi(Eigen::MatrixXd::Identity(2*robot.dimv(), 2*robot.dimv())),
      pi(Eigen::VectorXd::Zero(2*robot.dimv())),
      N(Eigen::MatrixXd::Zero(2*robot.dimv(), 2*robot.dimv())),
      n(Eigen::VectorXd::Zero(2*robot.dimv())) {
  }
    // : Pqq(Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv())),
    //   Pqv(Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv())),
    //   Pvq(Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv())),
    //   Pvv(Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv())),
    //   sq(Eigen::VectorXd::Zero(robot.dimv())),
    //   sv(Eigen::VectorXd::Zero(robot.dimv())),
    //   Pi(),
    //   pi(),
    //   N(),
    //   n() {
    // if (robot.max_point_contacts() > 0) {
    //   const int dimx = 2*robot.dimv();
    //   Pi.resize(dimx, dimx);
    //   Pi.setIdentity();
    //   pi.resize(dimx);
    //   pi.setZero();
    //   N.resize(dimx, dimx);
    //   N.setZero();
    //   n.resize(dimx);
    //   n.setZero();
    // }
  // }

  ///
  /// @brief Default constructor. 
  ///
  RiccatiFactorization()
    : Pqq(),
      Pqv(),
      Pvq(),
      Pvv(),
      sq(),
      sv(),
      Pi(),
      pi(),
      N(),
      n() {
  }

  ///
  /// @brief Destructor. 
  ///
  ~RiccatiFactorization() {
  }

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
  /// @brief Riccati factorization for pure-state equality constraints. 
  /// Size is 2 * Robot::dimv() x 2 * Robot::dimv() if the robot can have
  /// contacts. Size is 0 x 0 otherwise.
  ///
  Eigen::MatrixXd Pi;

  ///
  /// @brief Riccati factorization for pure-state equality constraints. Size is 
  /// 2 * Robot::dimv(). Size is 0 otherwise.
  ///
  Eigen::VectorXd pi;

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
  Eigen::VectorXd n;

  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

};

} // namespace idocp 

#endif // IDOCP_RICCATI_FACTORIZATION_HPP_ 