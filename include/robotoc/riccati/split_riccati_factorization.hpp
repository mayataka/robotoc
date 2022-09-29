#ifndef ROBOTOC_SPLIT_RICCATI_FACTORIZATION_HPP_
#define ROBOTOC_SPLIT_RICCATI_FACTORIZATION_HPP_

#include "Eigen/Core"

#include "robotoc/robot/robot.hpp"


namespace robotoc {

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
  SplitRiccatiFactorization(const Robot& robot);

  ///
  /// @brief Default constructor. 
  ///
  SplitRiccatiFactorization();

  ///
  /// @brief Destructor. 
  ///
  ~SplitRiccatiFactorization();

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
  /// @brief Riccati factorization vector w.r.t. the switching time. Size is 
  /// 2 * Robot::dimv().
  ///
  Eigen::VectorXd psi_x;

  ///
  /// @brief Riccati factorization vector w.r.t. the switching time. Size is 
  /// Robot::dimu().
  ///
  Eigen::VectorXd psi_u;

  ///
  /// @brief Riccati factorization vector w.r.t. the switching time. Size is 
  /// 2 * Robot::dimv().
  ///
  Eigen::VectorXd Psi;

  ///
  /// @brief Riccati factorization vector w.r.t. the switching time. Size is 
  /// 2 * Robot::dimv().
  ///
  Eigen::VectorXd phi_x;

  ///
  /// @brief Riccati factorization vector w.r.t. the switching time. Size is 
  /// Robot::dimu().
  ///
  Eigen::VectorXd phi_u;

  ///
  /// @brief Riccati factorization vector w.r.t. the switching time. Size is 
  /// 2 * Robot::dimv().
  ///
  Eigen::VectorXd Phi;

  ///
  /// @brief Riccati factorization w.r.t. the switching time. 
  ///
  double xi;

  ///
  /// @brief Riccati factorization w.r.t. the switching time. 
  ///
  double chi;

  ///
  /// @brief Riccati factorization w.r.t. the switching time. 
  ///
  double rho;

  ///
  /// @brief Riccati factorization w.r.t. the switching time. 
  ///
  double eta;

  ///
  /// @brief Riccati factorization w.r.t. the switching time. 
  ///
  double iota;

  void setConstraintDimension(const int dims=0);

  int dims() const;

  Eigen::Block<Eigen::MatrixXd> M();

  const Eigen::Block<const Eigen::MatrixXd> M() const;

  Eigen::VectorBlock<Eigen::VectorXd> m();

  const Eigen::VectorBlock<const Eigen::VectorXd> m() const;

  Eigen::VectorBlock<Eigen::VectorXd> mt();

  const Eigen::VectorBlock<const Eigen::VectorXd> mt() const;

    Eigen::VectorBlock<Eigen::VectorXd> mt_next();

  const Eigen::VectorBlock<const Eigen::VectorXd> mt_next() const;

  void setZero();

  void setRandom();

  ///
  /// @brief Checks the equivalence of two SplitRiccatiFactorization.
  /// @param[in] other object.
  /// @return true if this and other is same. false otherwise.
  ///
  bool isApprox(const SplitRiccatiFactorization& other) const;

  ///
  /// @brief Checks this object has at least one NaN.
  /// @return true if this has at least one NaN. false otherwise.
  ///
  bool hasNaN() const;

  static SplitRiccatiFactorization Random(const Robot& robot);

  ///
  /// @brief Displays the split Riccati factorization onto a ostream.
  ///
  void disp(std::ostream& os) const;

  friend std::ostream& operator<<(std::ostream& os, 
                                  const SplitRiccatiFactorization& riccati);

private:
  Eigen::MatrixXd M_full_;
  Eigen::VectorXd m_full_, mt_full_, mt_next_full_;
  int dimv_, dimx_, dims_;
};

} // namespace robotoc 

#include "robotoc/riccati/split_riccati_factorization.hxx"

#endif // ROBOTOC_SPLIT_RICCATI_FACTORIZATION_HPP_ 