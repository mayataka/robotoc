#ifndef ROBOTOC_SPLIT_KKT_RESIDUAL_HPP_ 
#define ROBOTOC_SPLIT_KKT_RESIDUAL_HPP_

#include <iostream>

#include "Eigen/Core"

#include "robotoc/robot/robot.hpp"
#include "robotoc/robot/contact_status.hpp"


namespace robotoc {

///
/// @class SplitKKTResidual
/// @brief KKT residual split into each time stage. 
///
class SplitKKTResidual {
public:
  ///
  /// @brief Construct a split KKT residual.
  /// @param[in] robot Robot model. 
  ///
  SplitKKTResidual(const Robot& robot);

  ///
  /// @brief Default constructor. 
  ///
  SplitKKTResidual();

  ///
  /// @brief Default destructor. 
  ///
  ~SplitKKTResidual() = default;

  ///
  /// @brief Default copy constructor. 
  ///
  SplitKKTResidual(const SplitKKTResidual&) = default;

  ///
  /// @brief Default copy operator. 
  ///
  SplitKKTResidual& operator=(const SplitKKTResidual&) = default;

  ///
  /// @brief Default move constructor. 
  ///
  SplitKKTResidual(SplitKKTResidual&&) noexcept = default;

  ///
  /// @brief Default move assign operator. 
  ///
  SplitKKTResidual& operator=(SplitKKTResidual&&) noexcept = default;

  ///
  /// @brief Sets contact status, i.e., set dimension of the contact forces.
  /// @param[in] dimf The dimension of the contact. Must be non-negative.
  ///
  void setContactDimension(const int dimf);

  ///
  /// @brief Sets the dimension of the switching constraint.
  /// @param[in] dims The dimension of the switching constraint. Must be non-negative.
  ///
  void setSwitchingConstraintDimension(const int dims);

  ///
  /// @brief Residual in the state equation. Size is 2 * Robot::dimv().
  ///
  Eigen::VectorXd Fx;

  ///
  /// @brief Residual in the state equation w.r.t. the configuration q.
  /// @return Reference to the residual in the state equation w.r.t. q. Size is 
  /// Robot::dimv().
  ///
  Eigen::VectorBlock<Eigen::VectorXd> Fq();

  ///
  /// @brief const version of SplitKKTResidual::Fq().
  ///
  const Eigen::VectorBlock<const Eigen::VectorXd> Fq() const;

  ///
  /// @brief Residual in the state equation w.r.t. the velocity v.
  /// @return Reference to the residual in the state equation w.r.t. v. Size is 
  /// Robot::dimv().
  ///
  Eigen::VectorBlock<Eigen::VectorXd> Fv();

  ///
  /// @brief const version of SplitKKTResidual::Fq().
  ///
  const Eigen::VectorBlock<const Eigen::VectorXd> Fv() const;

  ///
  /// @brief Residual in the switching constraint.
  /// @return Reference to the residual in the switching constraints. 
  /// Size is SplitKKTResidual::dims().
  ///
  Eigen::VectorBlock<Eigen::VectorXd> P();

  ///
  /// @brief const version of SplitKKTResidual::P().
  ///
  const Eigen::VectorBlock<const Eigen::VectorXd> P() const;

  ///
  /// @brief KKT Residual w.r.t. the state x. Size is 2 * Robot::dimv().
  ///
  Eigen::VectorXd lx;

  ///
  /// @brief KKT residual w.r.t. the configuration q. 
  /// @return Reference to the KKT residual w.r.t. q. Size is Robot::dimv().
  ///
  Eigen::VectorBlock<Eigen::VectorXd> lq();

  ///
  /// @brief const version of SplitKKTResidual::lq().
  ///
  const Eigen::VectorBlock<const Eigen::VectorXd> lq() const;

  ///
  /// @brief KKT residual w.r.t. the joint velocity v. 
  /// @return Reference to the KKT residual w.r.t. v. Size is Robot::dimv().
  ///
  Eigen::VectorBlock<Eigen::VectorXd> lv();

  ///
  /// @brief const version of SplitKKTResidual::lv().
  ///
  const Eigen::VectorBlock<const Eigen::VectorXd> lv() const;

  /// 
  /// @brief KKT residual w.r.t. the acceleration a. Size is Robot::dimv().
  /// 
  Eigen::VectorXd la;

  /// 
  /// @brief KKT residual w.r.t. the impact change in the velocity ddv. 
  /// Size is Robot::dimv().
  /// 
  Eigen::VectorXd ldv;

  /// 
  /// @brief KKT residual w.r.t. the control input torques u. Size is 
  /// Robot::dimu().
  /// 
  Eigen::VectorXd lu;

  ///
  /// @brief KKT residual w.r.t. the stack of the contact forces f. 
  /// @return Reference to the residual w.r.t. f. Size is 
  /// SplitKKTResidual::dimf().
  ///
  Eigen::VectorBlock<Eigen::VectorXd> lf();

  ///
  /// @brief const version of SplitKKTResidual::lf().
  ///
  const Eigen::VectorBlock<const Eigen::VectorXd> lf() const;

  ///
  /// @brief KKT residual w.r.t. the switching time, that is, this is the value
  /// of the Hamiltonian. 
  ///
  double h;

  ///
  /// @brief Returns the squared norm of the KKT residual, that is, 
  /// the primal and dual residual. 
  /// @return The squared norm of the KKT residual.
  ///
  double KKTError() const;

  ///
  /// @brief Returns the lp norm of the primal feasibility, i.e., the constraint 
  /// violation. Default norm is l1-norm. You can also specify l-infty norm by 
  /// passing Eigen::Infinity as the template parameter.
  /// @tparam p Index of norm. Default is 1 (l1-norm).
  /// @return The lp norm of the primal feasibility.
  ///
  template <int p=1>
  double primalFeasibility() const;

  ///
  /// @brief Returns the lp norm of the dual feasibility. Default norm is 
  /// l1-norm. You can also specify l-infty norm by passing Eigen::Infinity as 
  /// the template parameter.
  /// @tparam p Index of norm. Default is 1 (l1-norm).
  /// @return The lp norm of the dual feasibility.
  ///
  template <int p=1>
  double dualFeasibility() const;

  ///
  /// @brief Sets the split KKT residual zero.
  ///
  void setZero();

  ///
  /// @brief Returns the dimension of the stack of the contact forces at the 
  /// current contact status.
  /// @return Dimension of the stack of the contact forces.
  ///
  int dimf() const;

  ///
  /// @brief Returns the dimension of the stack of the contact forces at the 
  /// current contact status.
  /// @return Dimension of the stack of the contact forces.
  ///
  int dims() const;

  ///
  /// @brief Checks dimensional consistency of each component. 
  /// @return true if the dimension is consistent. false if not.
  ///
  bool isDimensionConsistent() const;

  ///
  /// @brief Checks the equivalence of two SplitKKTResidual.
  /// @param[in] other Other object.
  /// @return true if this and other is same. false otherwise.
  ///
  bool isApprox(const SplitKKTResidual& other) const;

  ///
  /// @brief Checks this has at least one NaN.
  /// @return true if this has at least one NaN. false otherwise.
  ///
  bool hasNaN() const;

  ///
  /// @brief Set by random value based on the current contact status. 
  ///
  void setRandom();

  ///
  /// @brief Set by random value. Contact status is reset.
  /// @param[in] contact_status Contact status.
  ///
  void setRandom(const ContactStatus& contact_status);

  ///
  /// @brief Set by random value. Contact status is reset.
  /// @param[in] impact_status Contact status.
  ///
  void setRandom(const ImpactStatus& impact_status);

  ///
  /// @brief Generates split KKT residual filled randomly.
  /// @return Split KKT residual filled randomly.
  /// @param[in] robot Robot model. 
  ///
  static SplitKKTResidual Random(const Robot& robot);

  ///
  /// @brief Generates split KKT residual filled randomly.
  /// @return Split KKT residual filled randomly.
  /// @param[in] robot Robot model. 
  /// @param[in] contact_status Contact status.
  ///
  static SplitKKTResidual Random(const Robot& robot, 
                                 const ContactStatus& contact_status);

  ///
  /// @brief Generates split KKT residual filled randomly.
  /// @return Split KKT residual filled randomly.
  /// @param[in] robot Robot model. 
  /// @param[in] impact_status Contact status.
  ///
  static SplitKKTResidual Random(const Robot& robot, 
                                 const ImpactStatus& impact_status);

  ///
  /// @brief Displays the split KKT residual onto a ostream.
  ///
  void disp(std::ostream& os) const;

  friend std::ostream& operator<<(std::ostream& os, 
                                  const SplitKKTResidual& kkt_residual);

private:
  Eigen::VectorXd P_full_, lf_full_;
  int dimv_, dimu_, dimf_, dims_;

};

} // namespace robotoc 

#include "robotoc/core/split_kkt_residual.hxx"

#endif // ROBOTOC_SPLIT_KKT_RESIDUAL_HPP_