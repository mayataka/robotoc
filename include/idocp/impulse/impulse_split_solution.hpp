#ifndef IDOCP_IMPULSE_SPLIT_SOLUTION_HPP_ 
#define IDOCP_IMPULSE_SPLIT_SOLUTION_HPP_

#include <vector>

#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/robot/impulse_status.hpp"
#include "idocp/impulse/impulse_split_direction.hpp"


namespace idocp {

///
/// @class ImpulseSplitSolution
/// @brief Solution of the optimal control problem at the impulse stage. 
///
class ImpulseSplitSolution {
public:
  ///
  /// @brief Construct a split solution.
  /// @param[in] robot Robot model. Must be initialized by URDF or XML.
  ///
  ImpulseSplitSolution(const Robot& robot);

  ///
  /// @brief Default constructor. 
  ///
  ImpulseSplitSolution();

  ///
  /// @brief Destructor. 
  ///
  ~ImpulseSplitSolution();

  ///
  /// @brief Default copy constructor. 
  ///
  ImpulseSplitSolution(const ImpulseSplitSolution&) = default;

  ///
  /// @brief Default copy assign operator. 
  ///
  ImpulseSplitSolution& operator=(const ImpulseSplitSolution&) = default;

  ///
  /// @brief Default move constructor. 
  ///
  ImpulseSplitSolution(ImpulseSplitSolution&&) noexcept = default;

  ///
  /// @brief Default move assign operator. 
  ///
  ImpulseSplitSolution& operator=(ImpulseSplitSolution&&) noexcept = default;

  ///
  /// @brief Set impulse status, i.e., set dimension of the impulse.
  /// @param[in] impulse_status Impulse status.
  ///
  void setImpulseStatus(const ImpulseStatus& impulse_status);

  ///
  /// @brief Stack of active impulse forces. Size is ImpulseStatus::dimp().
  /// @return Reference to the stack of active impulse forces.
  ///
  Eigen::VectorBlock<Eigen::VectorXd> f_stack();

  ///
  /// @brief const version of ImpulseSplitSolution::f_stack().
  ///
  const Eigen::VectorBlock<const Eigen::VectorXd> f_stack() const;

  ///
  /// @brief Set ImpulseSplitSolution::f_stack() from ImpulseSplitSolution::f.
  ///
  void set_f_stack();

  ///
  /// @brief Set ImpulseSplitSolution::f from ImpulseSplitSolution::f_stack().
  ///
  void set_f_vector();

  ///
  /// @brief Stack of Lagrange multiplier with respect to impulse velocity 
  /// constraints that is active at the current impulse status. Size is 
  /// ImpulseSplitSolution::dimf().
  /// @return Reference to the stack of Lagrange multiplier with respect to 
  /// impulse velocity constraints.
  ///
  Eigen::VectorBlock<Eigen::VectorXd> mu_stack();

  ///
  /// @brief const version of ImpulseSplitSolution::mu_stack().
  ///
  const Eigen::VectorBlock<const Eigen::VectorXd> mu_stack() const;

  ///
  /// @brief Set ImpulseSplitSolution::mu_stack() from ImpulseSplitSolution::mu. 
  ///
  void set_mu_stack();

  ///
  /// @brief Set ImpulseSplitSolution::mu from ImpulseSplitSolution::mu_stack(). 
  ///
  void set_mu_vector();

  ///
  /// @brief Stack of Lagrange multiplier with respect to impulse position 
  /// constraints that is active at the current impulse status. Size is 
  /// ImpulseSplitSolution::dimf().
  /// @return Reference to the stack of Lagrange multiplier with respect to 
  /// impulse position constraints.
  ///
  Eigen::VectorBlock<Eigen::VectorXd> xi_stack();

  ///
  /// @brief const version of ImpulseSplitSolution::xi_stack().
  ///
  const Eigen::VectorBlock<const Eigen::VectorXd> xi_stack() const;

  ///
  /// @brief Set ImpulseSplitSolution::xi_stack() from ImpulseSplitSolution::mu. 
  ///
  void set_xi_stack();

  ///
  /// @brief Set ImpulseSplitSolution::xi from ImpulseSplitSolution::mu_stack(). 
  ///
  void set_xi_vector();

  ///
  /// @brief Returns the dimension of the stack of impulse forces at the current 
  /// impulse status.
  /// @return Dimension of impulse forces.
  ///
  int dimf() const;

  ///
  /// @brief Lagrange multiplier with respect to the transition of 
  /// ImpulseSplitSolution::q. Size is ImpulseStatus::dimp().
  ///
  Eigen::VectorXd lmd;

  ///
  /// @brief Lagrange multiplier with respect to the transition of  
  /// ImpulseSplitSolution::v. Size is Robot::dimv().
  ///
  Eigen::VectorXd gmm;

  ///
  /// @brief Configuration. Size is Robot::dimq().
  ///
  Eigen::VectorXd q;

  ///
  /// @brief Generalized velocity. Size is Robot::dimv().
  ///
  Eigen::VectorXd v;

  ///
  /// @brief Impulse change in the generalized velocity. 
  /// Size is Robot::dimv().
  ///
  Eigen::VectorXd dv;

  ///
  /// @brief Impulse forces. 
  /// Size is Robot::max_point_contacts().
  ///
  std::vector<Eigen::Vector3d> f;

  ///
  /// @brief Lagrange multiplier with respect to impulse dynamics. Size is 
  /// Robot::dimv().
  ///
  Eigen::VectorXd beta;

  ///
  /// @brief Lagrange multiplier with respect to impulse velocity constraint. 
  /// Size is Robot::max_point_contacts().
  ///
  std::vector<Eigen::Vector3d> mu;

  ///
  /// @brief Lagrange multiplier with respect to impulse position constraint. 
  /// Size is Robot::max_point_contacts().
  ///
  std::vector<Eigen::Vector3d> xi;

  ///
  /// @brief Return true if a Impulse is active and false if not.
  /// @param[in] contact_index Index of a contact at impulse. 
  /// @return true if a Impulse is active and false if not. 
  ///
  bool isImpulseActive(const int contact_index) const;

  ///
  /// @brief Integrates the solution based on step size and direction. 
  /// @param[in] robot Robot model.
  /// @param[in] step_size Step size.
  /// @param[in] d Split direction.
  ///
  void integrate(const Robot& robot, const double step_size, 
                 const ImpulseSplitDirection& d);

  ///
  /// @brief Return true if two ImpulseSplitSolution have the same value and  
  /// false if not. 
  /// @param[in] other Impulse split solution that is compared with this object.
  ///
  bool isApprox(const ImpulseSplitSolution& other) const;

  ///
  /// @brief Set each component vector by random value based on the current 
  /// impulse status. 
  /// @param[in] robot Robot model.
  ///
  void setRandom(const Robot& robot);

  ///
  /// @brief Set each component vector by random value. Impulse status is reset.
  /// @param[in] robot Robot model.
  /// @param[in] impulse_status Impulse status.
  ///
  void setRandom(const Robot& robot, const ImpulseStatus& impulse_status);

  ///
  /// @brief Generates impulse split solution filled randomly.
  /// @return Impulse split solution filled randomly.
  /// @param[in] robot Robot model. 
  ///
  static ImpulseSplitSolution Random(const Robot& robot);

  ///
  /// @brief Generates impulse split solution filled randomly.
  /// @return Impulse split solution filled randomly.
  /// @param[in] robot Robot model. 
  /// @param[in] impulse_status Impulse status.
  ///
  static ImpulseSplitSolution Random(const Robot& robot, 
                                     const ImpulseStatus& impulse_status);

  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

private:
  Eigen::VectorXd mu_stack_, f_stack_, xi_stack_;
  bool has_floating_base_;
  std::vector<bool> is_impulse_active_;
  int dimf_;

};

} // namespace idocp 

#include "idocp/impulse/impulse_split_solution.hxx"

#endif // IDOCP_IMPULSE_SPLIT_SOLUTION_HPP_ 