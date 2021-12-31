#ifndef ROBOTOC_IMPULSE_WRENCH_FRICTION_CONE_HPP_
#define ROBOTOC_IMPULSE_WRENCH_FRICTION_CONE_HPP_

#include "Eigen/Core"

#include "robotoc/robot/robot.hpp"
#include "robotoc/robot/impulse_status.hpp"
#include "robotoc/impulse/impulse_split_solution.hpp"
#include "robotoc/impulse/impulse_split_direction.hpp"
#include "robotoc/constraints/impulse_constraint_component_base.hpp"
#include "robotoc/constraints/constraint_component_data.hpp"
#include "robotoc/impulse/impulse_split_kkt_residual.hpp"
#include "robotoc/impulse/impulse_split_kkt_matrix.hpp"


namespace robotoc {

///
/// @class ImpulseWrenchFrictionCone
/// @brief Constraint on the wrench firction cone for surface contacts.
///
class ImpulseWrenchFrictionCone final : public ImpulseConstraintComponentBase {
public:
  ///
  /// @brief Constructor. 
  /// @param[in] robot Robot model.
  /// @param[in] mu Friction coefficient. Must be positive.
  /// @param[in] X A length of the rectangular. Must be positive.
  /// @param[in] Y A length of the rectangular. Must be positive.
  ///
  ImpulseWrenchFrictionCone(const Robot& robot, const double mu,
                            const double X, const double Y);

  ///
  /// @brief Default constructor. 
  ///
  ImpulseWrenchFrictionCone();

  ///
  /// @brief Destructor. 
  ///
  ~ImpulseWrenchFrictionCone();

  ///
  /// @brief Default copy constructor. 
  ///
  ImpulseWrenchFrictionCone(const ImpulseWrenchFrictionCone&) = default;

  ///
  /// @brief Default copy operator. 
  ///
  ImpulseWrenchFrictionCone& operator=(const ImpulseWrenchFrictionCone&) = default;

  ///
  /// @brief Default move constructor. 
  ///
  ImpulseWrenchFrictionCone(ImpulseWrenchFrictionCone&&) noexcept = default;

  ///
  /// @brief Default move assign operator. 
  ///
  ImpulseWrenchFrictionCone& operator=(ImpulseWrenchFrictionCone&&) noexcept = default;

  ///
  /// @brief Sets the friction coefficient. 
  /// @param[in] mu Friction coefficient. Must be positive.
  ///
  void setFrictionCoefficient(const double mu);

  ///
  /// @param[in] X A length of the rectangular. Must be positive.
  /// @param[in] Y A length of the rectangular. Must be positive.
  ///
  void setRectangular(const double X, const double Y);

  KinematicsLevel kinematicsLevel() const override;

  void allocateExtraData(ConstraintComponentData& data) const override;

  bool isFeasible(Robot& robot, const ImpulseStatus& impulse_status, 
                  ConstraintComponentData& data, 
                  const ImpulseSplitSolution& s) const override;

  void setSlack(Robot& robot, const ImpulseStatus& impulse_status, 
                ConstraintComponentData& data, 
                const ImpulseSplitSolution& s) const override;

  void evalConstraint(Robot& robot, const ImpulseStatus& impulse_status, 
                      ConstraintComponentData& data, 
                      const ImpulseSplitSolution& s) const override;

  void evalDerivatives(Robot& robot, const ImpulseStatus& impulse_status, 
                       ConstraintComponentData& data, 
                       const ImpulseSplitSolution& s,
                       ImpulseSplitKKTResidual& kkt_residual) const override;

  void condenseSlackAndDual(const ImpulseStatus& impulse_status,
                            ConstraintComponentData& data, 
                            ImpulseSplitKKTMatrix& kkt_matrix,
                            ImpulseSplitKKTResidual& kkt_residual) const override;

  void expandSlackAndDual(const ImpulseStatus& impulse_status, 
                          ConstraintComponentData& data, 
                          const ImpulseSplitDirection& d) const override; 

  int dimc() const override;

  EIGEN_MAKE_ALIGNED_OPERATOR_NEW 

private:
  int dimv_, dimc_, max_num_contacts_;
  std::vector<int> contact_frame_;
  std::vector<ContactType> contact_types_;
  double mu_, X_, Y_;
  Eigen::MatrixXd cone_;

  void setCone(const double mu, const double X, const double Y);

};

} // namespace robotoc

#endif // ROBOTOC_IMPULSE_WRENCH_FRICTION_CONE_HPP_ 