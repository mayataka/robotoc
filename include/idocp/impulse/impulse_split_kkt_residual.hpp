#ifndef IDOCP_IMPULSE_SPLIT_KKT_RESIDUAL_HPP_ 
#define IDOCP_IMPULSE_SPLIT_KKT_RESIDUAL_HPP_

#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/robot/impulse_status.hpp"


namespace idocp {

///
/// @class ImpulseSplitKKTResidual
/// @brief KKT residual split into an impulse stage. 
///
class ImpulseSplitKKTResidual {
public:
  ///
  /// @brief Construct a split KKT residual.
  /// @param[in] robot Robot model. 
  ///
  ImpulseSplitKKTResidual(const Robot& robot);

  ///
  /// @brief Default constructor. 
  ///
  ImpulseSplitKKTResidual();

  ///
  /// @brief Destructor. 
  ///
  ~ImpulseSplitKKTResidual();

  ///
  /// @brief Default copy constructor. 
  ///
  ImpulseSplitKKTResidual(const ImpulseSplitKKTResidual&) = default;

  ///
  /// @brief Default copy operator. 
  ///
  ImpulseSplitKKTResidual& operator=(const ImpulseSplitKKTResidual&) = default;

  ///
  /// @brief Default move constructor. 
  ///
  ImpulseSplitKKTResidual(ImpulseSplitKKTResidual&&) noexcept = default;

  ///
  /// @brief Default move assign operator. 
  ///
  ImpulseSplitKKTResidual& operator=(
      ImpulseSplitKKTResidual&&) noexcept = default;

  ///
  /// @brief Sets the impulse status, i.e., set the dimension of the impulse.
  /// @param[in] impulse_status Impulse status.
  ///
  void setImpulseStatus(const ImpulseStatus& impulse_status);

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
  /// @brief const version of ImpulseSplitKKTResidual::Fq().
  ///
  const Eigen::VectorBlock<const Eigen::VectorXd> Fq() const;

  ///
  /// @brief Residual in the state equation w.r.t. the joint velocity v.
  /// @return Reference to the residual in the state equation w.r.t. v. Size is 
  /// Robot::dimv().
  ///
  Eigen::VectorBlock<Eigen::VectorXd> Fv();

  ///
  /// @brief const version of ImpulseSplitKKTResidual::Fv().
  ///
  const Eigen::VectorBlock<const Eigen::VectorXd> Fv() const;

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
  /// @brief const version of ImpulseSplitKKTResidual::lq().
  ///
  const Eigen::VectorBlock<const Eigen::VectorXd> lq() const;

  ///
  /// @brief KKT residual w.r.t. the joint velocity v. 
  /// @return Reference to the KKT residual w.r.t. v. Size is Robot::dimv().
  ///
  Eigen::VectorBlock<Eigen::VectorXd> lv();

  ///
  /// @brief const version of ImpulseSplitKKTResidual::lv().
  ///
  const Eigen::VectorBlock<const Eigen::VectorXd> lv() const;

  /// 
  /// @brief KKT residual w.r.t. the impulse change in the velocity ddv. 
  /// Size is Robot::dimv().
  /// 
  Eigen::VectorXd ldv;

  ///
  /// @brief KKT residual w.r.t. the stack of the impulse forces f. 
  /// @return Reference to the residual w.r.t. f. Size is 
  /// ImpulseSplitKKTResidual::dimi().
  ///
  Eigen::VectorBlock<Eigen::VectorXd> lf();

  ///
  /// @brief const version of ImpulseSplitKKTResidual::lf().
  ///
  const Eigen::VectorBlock<const Eigen::VectorXd> lf() const;

  ///
  /// @brief Returns the squared norm of the KKT residual, that is, 
  /// the primal and dual residual. 
  /// @return The squared norm of the KKT residual.
  ///
  double KKTError() const;

  ///
  /// @brief Returns the l1 norm of the constraint violation, that is,
  /// the primal residual in the state equation. 
  /// @return The l1 norm of the constraint violation.
  ///
  double constraintViolation() const;

  ///
  /// @brief Sets the split KKT residual zero.
  ///
  void setZero();

  ///
  /// @brief Returns the dimension of the stack of the impluse forces at the 
  /// current impulse status.
  /// @return Dimension of the stack of the impulse forces.
  ///
  int dimi() const;

  ///
  /// @brief Checks dimensional consistency of each component. 
  /// @return true if the dimension is consistent. false if not.
  ///
  bool isDimensionConsistent() const;

  ///
  /// @brief Checks the equivalence of two ImpulseSplitKKTResidual.
  /// @param[in] other Other object.
  /// @return true if this and other is same. false otherwise.
  ///
  bool isApprox(const ImpulseSplitKKTResidual& other) const;

  ///
  /// @brief Checks if this has at least one NaN.
  /// @return true if this has at least one NaN. false otherwise.
  ///
  bool hasNaN() const;

  ///
  /// @brief Set by random value based on the current impulse status. 
  ///
  void setRandom();

  ///
  /// @brief Set by random value. Impulse status is reset.
  /// @param[in] impulse_status Impulse status.
  ///
  void setRandom(const ImpulseStatus& impulse_status);

  ///
  /// @brief Generates impulse split KKT residual filled randomly.
  /// @return Impulse split KKT residual filled randomly.
  /// @param[in] robot Robot model. 
  ///
  static ImpulseSplitKKTResidual Random(const Robot& robot);

  ///
  /// @brief Generates impulse split KKT residual filled randomly.
  /// @return Impulse split KKT residual filled randomly.
  /// @param[in] robot Robot model. 
  /// @param[in] impulse_status Impulse status.
  ///
  static ImpulseSplitKKTResidual Random(const Robot& robot, 
                                        const ImpulseStatus& impulse_status);

private:
  Eigen::VectorXd lf_full_;
  int dimv_, dimi_;

};

} // namespace idocp 

#include "idocp/impulse/impulse_split_kkt_residual.hxx"

#endif // IDOCP_IMPULSE_SPLIT_KKT_RESIDUAL_HPP_ 