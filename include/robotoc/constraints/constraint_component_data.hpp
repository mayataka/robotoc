#ifndef ROBOTOC_CONSTRAINT_COMPONENT_DATA_HPP_
#define ROBOTOC_CONSTRAINT_COMPONENT_DATA_HPP_

#include <vector>

#include "Eigen/Core"


namespace robotoc {

///
/// @class ConstraintComponentData
/// @brief Data used in constraint components. Composed by slack, 
/// dual (Lagrange multiplier), primal residual, complementary slackness between 
/// the slack and dual, and directions of slack and dual.
///
class ConstraintComponentData {
public:
  ///
  /// @brief Constructor. 
  /// @param[in] dimc Dimension of the constraint component. Must be positive.
  /// @param[in] barrier_param Barrier parameter. Must be positive. Should be small.
  /// Only used to initialize the slack and dual variables.
  ///
  ConstraintComponentData(const int dimc, const double barrier_param);

  ///
  /// @brief Default constructor. 
  ///
  ConstraintComponentData();

  ///
  /// @brief Default destructor. 
  ///
  ~ConstraintComponentData() = default;

  ///
  /// @brief Default copy constructor. 
  ///
  ConstraintComponentData(const ConstraintComponentData&) = default;

  ///
  /// @brief Default copy operator. 
  ///
  ConstraintComponentData& operator=(const ConstraintComponentData&) = default;

  ///
  /// @brief Default move constructor. 
  ///
  ConstraintComponentData(ConstraintComponentData&&) noexcept = default;

  ///
  /// @brief Default move assign operator. 
  ///
  ConstraintComponentData& operator=(ConstraintComponentData&&) noexcept = default;

  ///
  /// @brief Slack variable of the constraint. Size is 
  /// ConstraintComponentData::dimc(). All elements must be positive.
  ///
  Eigen::VectorXd slack;

  ///
  /// @brief Dual variable (Lagrange multiplier) of the constraint. Size is 
  /// ConstraintComponentData::dimc(). All elements must be positive.
  ///
  Eigen::VectorXd dual;

  ///
  /// @brief Primal residual of the constraint. Size is 
  /// ConstraintComponentData::dimc(). 
  ///
  Eigen::VectorXd residual;

  ///
  /// @brief Residual in the complementary slackness between slack and dual. 
  /// Size is ConstraintComponentData::dimc(). 
  ///
  Eigen::VectorXd cmpl;

  ///
  /// @brief Newton direction of the slack. Size is 
  /// ConstraintComponentData::dimc(). 
  ///
  Eigen::VectorXd dslack;

  ///
  /// @brief Newton direction of the dual. Size is 
  /// ConstraintComponentData::dimc(). 
  ///
  Eigen::VectorXd ddual;

  ///
  /// @brief Used in condensing of slack and dual. Size is 
  /// ConstraintComponentData::dimc(). 
  ///
  Eigen::VectorXd cond;

  ///
  /// @brief Value of the log berrier function of the slack variable.
  ///
  double log_barrier;

  ///
  /// @brief std vector of Eigen::VectorXd used to store residual temporaly. 
  /// Only be allocated in ConstraintComponentBase::allocateExtraData().
  ///
  std::vector<Eigen::VectorXd> r;

  ///
  /// @brief std vector of Eigen::MatrixXd used to store Jacobian temporaly. 
  /// Only be allocated in ConstraintComponentBase::allocateExtraData().
  ///
  std::vector<Eigen::MatrixXd> J;

  ///
  /// @brief Copies the slack and dual variables from another constraint 
  /// component data. this->dimc() and other.dimc() must be the same.
  /// @param[in] other Another constraint component data. 
  ///
  void copySlackAndDual(const ConstraintComponentData& other);

  ///
  /// @brief Returns the squared norm of the KKT reisdual, that is, the sum of
  /// the squared norm of the primal residual and complementary slackness of 
  /// the constraint. 
  /// @return Squared norm of the KKT residual. 
  ///
  double KKTError() const {
    return (residual.squaredNorm() + cmpl.squaredNorm());
  }

  ///
  /// @brief Returns the lp norm of the primal feasibility, i.e., the constraint 
  /// violation. Default norm is l1-norm. You can also specify l-infty norm by 
  /// passing Eigen::Infinity as the template parameter.
  /// @tparam p Index of norm. Default is 1 (l1-norm).
  /// @return The lp norm of the primal feasibility.
  ///
  template <int p=1>
  double primalFeasibility() const {
    return residual.template lpNorm<p>();
  }

  ///
  /// @brief Returns the lp norm of the dual feasibility. Default norm is 
  /// l1-norm. You can also specify l-infty norm by passing Eigen::Infinity as 
  /// the template parameter.
  /// @tparam p Index of norm. Default is 1 (l1-norm).
  /// @return The lp norm of the dual feasibility.
  ///
  template <int p=1>
  double dualFeasibility() const {
    return cmpl.template lpNorm<p>();
  }

  ///
  /// @brief Resizes the constraint. 
  /// @param[in] dimc The new size. 
  ///
  void resize(const int dimc);

  ///
  /// @brief Dimension of the constraint. 
  /// @return Dimension of the constraint. 
  ///
  int dimc() const {
    return dimc_;
  }

  ///
  /// @brief Check whether dimensions of slack, dual, residual, cmpl, 
  /// dslack, ddual are ConstraintComponentData::dimc(). 
  /// @return Dimension of the constraint. 
  ///
  bool checkDimensionalConsistency() const;

  ///
  /// @brief Checks the equivalence of two ConstraintComponentData.
  /// @param[in] other Other object.
  /// @return true if this and other is same. false otherwise.
  ///
  bool isApprox(const ConstraintComponentData& other) const;

private:
  int dimc_;

};

} // namespace robotoc

#endif // ROBOTOC_CONSTRAINT_COMPONENT_DATA_HPP_