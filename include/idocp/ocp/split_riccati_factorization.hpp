#ifndef IDOCP_SPLIT_RICCATI_FACTORIZATION_HPP_
#define IDOCP_SPLIT_RICCATI_FACTORIZATION_HPP_

#include "Eigen/Core"

#include "idocp/robot/robot.hpp"


namespace idocp {

///
/// @class SplitRiccatiFactorization
/// @brief Riccati factorized matrix and vector split into a time stage.
///
struct SplitRiccatiFactorization {
public:
  ///
  /// @brief Allocate Riccati factorization matrix and vector.
  /// @param[in] robot Robot model. Must be initialized by URDF or XML.
  ///
  SplitRiccatiFactorization(const Robot& robot)
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
    // if (robot.maxPointContacts() > 0) {
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
  SplitRiccatiFactorization()
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
  ~SplitRiccatiFactorization() {
  }

  ///
  /// @brief Default copy constructor. 
  ///
  SplitRiccatiFactorization(const SplitRiccatiFactorization&) = default;

  ///
  /// @brief Default copy operator. 
  ///
  SplitRiccatiFactorization& operator=(
      const SplitRiccatiFactorization&) = default;

  ///
  /// @brief Default move constructor. 
  ///
  SplitRiccatiFactorization(
      SplitRiccatiFactorization&&) noexcept = default;

  ///
  /// @brief Default move assign operator. 
  ///
  SplitRiccatiFactorization& operator=(
      SplitRiccatiFactorization&&) noexcept = default;

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

  ///
  /// @brief Chech the equivalence of two SplitRiccatiFactorization.
  /// @param[in] other object.
  /// @return true if this and other is same. false otherwise.
  ///
  bool isApprox(const SplitRiccatiFactorization& other) const {
    if (!Pqq.isApprox(other.Pqq)) return false;
    if (!Pqv.isApprox(other.Pqv)) return false;
    if (!Pvq.isApprox(other.Pvq)) return false;
    if (!Pvv.isApprox(other.Pvv)) return false;
    if (!sq.isApprox(other.sq)) return false;
    if (!sv.isApprox(other.sv)) return false;
    if (!Pi.isApprox(other.Pi)) return false;
    if (!pi.isApprox(other.pi)) return false;
    if (!N.isApprox(other.N)) return false;
    if (!n.isApprox(other.n)) return false;
    return true;
  }

  ///
  /// @brief Chech this has at least one NaN.
  /// @return true if this has at least one NaN. false otherwise.
  ///
  bool hasNaN() const {
    if (Pqq.hasNaN()) return true;
    if (Pqv.hasNaN()) return true;
    if (Pvq.hasNaN()) return true;
    if (Pvv.hasNaN()) return true;
    if (sq.hasNaN()) return true;
    if (sv.hasNaN()) return true;
    if (Pi.hasNaN()) return true;
    if (pi.hasNaN()) return true;
    if (N.hasNaN()) return true;
    if (n.hasNaN()) return true;
    return false;
  }

  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

};

} // namespace idocp 

#endif // IDOCP_SPLIT_RICCATI_FACTORIZATION_HPP_ 