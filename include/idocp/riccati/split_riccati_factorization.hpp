#ifndef IDOCP_SPLIT_RICCATI_FACTORIZATION_HPP_
#define IDOCP_SPLIT_RICCATI_FACTORIZATION_HPP_

#include "Eigen/Core"

#include "idocp/robot/robot.hpp"


namespace idocp {

///
/// @class SplitRiccatiFactorization
/// @brief Riccati factorization matrix and vector for a time stage.
///
struct SplitRiccatiFactorization {
public:
  ///
  /// @brief Constructs Riccati factorization matrix and vector.
  /// @param[in] robot Robot model. 
  ///
  SplitRiccatiFactorization(const Robot& robot)
    : Pqq(Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv())),
      Pqv(Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv())),
      Pvq(Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv())),
      Pvv(Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv())),
      sq(Eigen::VectorXd::Zero(robot.dimv())),
      sv(Eigen::VectorXd::Zero(robot.dimv())) {
  }

  ///
  /// @brief Default constructor. 
  ///
  SplitRiccatiFactorization()
    : Pqq(),
      Pqv(),
      Pvq(),
      Pvv(),
      sq(),
      sv() {
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
  /// @brief Checks the equivalence of two SplitRiccatiFactorization.
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
    return true;
  }

  ///
  /// @brief Checks this object has at least one NaN.
  /// @return true if this has at least one NaN. false otherwise.
  ///
  bool hasNaN() const {
    if (Pqq.hasNaN()) return true;
    if (Pqv.hasNaN()) return true;
    if (Pvq.hasNaN()) return true;
    if (Pvv.hasNaN()) return true;
    if (sq.hasNaN()) return true;
    if (sv.hasNaN()) return true;
    return false;
  }

};

} // namespace idocp 

#endif // IDOCP_SPLIT_RICCATI_FACTORIZATION_HPP_ 