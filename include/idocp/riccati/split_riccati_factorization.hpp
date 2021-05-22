#ifndef IDOCP_SPLIT_RICCATI_FACTORIZATION_HPP_
#define IDOCP_SPLIT_RICCATI_FACTORIZATION_HPP_

#include "Eigen/Core"

#include "idocp/robot/robot.hpp"


namespace idocp {

///
/// @class SplitRiccatiFactorization
/// @brief Riccati factorization matrix and vector for a time stage.
///
class SplitRiccatiFactorization {
public:
  ///
  /// @brief Constructs Riccati factorization matrix and vector.
  /// @param[in] robot Robot model. 
  ///
  SplitRiccatiFactorization(const Robot& robot)
    : P(Eigen::MatrixXd::Zero(2*robot.dimv(), 2*robot.dimv())),
      s(Eigen::VectorXd::Zero(2*robot.dimv())),
      dimv_(robot.dimv()),
      dimx_(2*robot.dimv()) {
  }

  ///
  /// @brief Default constructor. 
  ///
  SplitRiccatiFactorization()
    : P(),
      s(),
      dimv_(0),
      dimx_(0) {
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
  /// 2 * Robot::dimv() x 2 * Robot::dimv().
  ///
  Eigen::MatrixXd P;

  ///
  /// @brief Riccati factorization vector. Size is 2 * Robot::dimv().
  ///
  Eigen::VectorXd s;

  Eigen::Block<Eigen::MatrixXd> Pqq() {
    return P.topLeftCorner(dimv_, dimv_); 
  }

  const Eigen::Block<const Eigen::MatrixXd> Pqq() const {
    return P.topLeftCorner(dimv_, dimv_); 
  }

  Eigen::Block<Eigen::MatrixXd> Pqv() {
    return P.topRightCorner(dimv_, dimv_); 
  }

  const Eigen::Block<const Eigen::MatrixXd> Pqv() const {
    return P.topRightCorner(dimv_, dimv_); 
  }

  Eigen::Block<Eigen::MatrixXd> Pvq() {
    return P.bottomLeftCorner(dimv_, dimv_); 
  }

  const Eigen::Block<const Eigen::MatrixXd> Pvq() const {
    return P.bottomLeftCorner(dimv_, dimv_); 
  }

  Eigen::Block<Eigen::MatrixXd> Pvv() {
    return P.bottomRightCorner(dimv_, dimv_); 
  }

  const Eigen::Block<const Eigen::MatrixXd> Pvv() const {
    return P.bottomRightCorner(dimv_, dimv_); 
  }

  Eigen::VectorBlock<Eigen::VectorXd> sq() {
    return s.head(dimv_);
  }

  const Eigen::VectorBlock<const Eigen::VectorXd> sq() const {
    return s.head(dimv_);
  }

  Eigen::VectorBlock<Eigen::VectorXd> sv() {
    return s.tail(dimv_);
  }

  const Eigen::VectorBlock<const Eigen::VectorXd> sv() const {
    return s.tail(dimv_);
  }

  ///
  /// @brief Checks the equivalence of two SplitRiccatiFactorization.
  /// @param[in] other object.
  /// @return true if this and other is same. false otherwise.
  ///
  bool isApprox(const SplitRiccatiFactorization& other) const {
    if (!P.isApprox(other.P)) return false;
    if (!s.isApprox(other.s)) return false;
    return true;
  }

  ///
  /// @brief Checks this object has at least one NaN.
  /// @return true if this has at least one NaN. false otherwise.
  ///
  bool hasNaN() const {
    if (P.hasNaN()) return true;
    if (s.hasNaN()) return true;
    return false;
  }

private:
  int dimv_, dimx_;

};

} // namespace idocp 

#endif // IDOCP_SPLIT_RICCATI_FACTORIZATION_HPP_ 