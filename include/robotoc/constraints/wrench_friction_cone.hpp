#ifndef ROBOTOC_WRENCH_FRICTION_CONE_HPP_
#define ROBOTOC_WRENCH_FRICTION_CONE_HPP_

#include "Eigen/Core"

#include "robotoc/robot/robot.hpp"
#include "robotoc/robot/contact_status.hpp"
#include "robotoc/ocp/split_solution.hpp"
#include "robotoc/ocp/split_direction.hpp"
#include "robotoc/constraints/constraint_component_base.hpp"
#include "robotoc/constraints/constraint_component_data.hpp"
#include "robotoc/ocp/split_kkt_residual.hpp"
#include "robotoc/ocp/split_kkt_matrix.hpp"


namespace robotoc {

///
/// @class WrenchFrictionCone
/// @brief Constraint on the wrench firction cone for surface contacts.
///
class WrenchFrictionCone final : public ConstraintComponentBase {
public:
  ///
  /// @brief Constructor. 
  /// @param[in] robot Robot model.
  /// @param[in] mu Friction coefficient. Must be positive.
  /// @param[in] X A length of the rectangular. Must be positive.
  /// @param[in] Y A length of the rectangular. Must be positive.
  ///
  WrenchFrictionCone(const Robot& robot, const double mu,
                     const double X, const double Y);

  ///
  /// @brief Default constructor. 
  ///
  WrenchFrictionCone();

  ///
  /// @brief Destructor. 
  ///
  ~WrenchFrictionCone();

  ///
  /// @brief Default copy constructor. 
  ///
  WrenchFrictionCone(const WrenchFrictionCone&) = default;

  ///
  /// @brief Default copy operator. 
  ///
  WrenchFrictionCone& operator=(const WrenchFrictionCone&) = default;

  ///
  /// @brief Default move constructor. 
  ///
  WrenchFrictionCone(WrenchFrictionCone&&) noexcept = default;

  ///
  /// @brief Default move assign operator. 
  ///
  WrenchFrictionCone& operator=(WrenchFrictionCone&&) noexcept = default;

  DEFINE_DEFAULT_CLONE_CONSTRAINT_COMPONENT(WrenchFrictionCone)

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

  bool useKinematics() const override;

  KinematicsLevel kinematicsLevel() const override;

  void allocateExtraData(ConstraintComponentData& data) const override;

  bool isFeasible(Robot& robot, const ContactStatus& contact_status, 
                  ConstraintComponentData& data, 
                  const SplitSolution& s) const override;

  void setSlack(Robot& robot, const ContactStatus& contact_status, 
                ConstraintComponentData& data, 
                const SplitSolution& s) const override;

  void evalConstraint(Robot& robot, const ContactStatus& contact_status, 
                      ConstraintComponentData& data, 
                      const SplitSolution& s) const override;

  void evalDerivatives(Robot& robot, const ContactStatus& contact_status, 
                       ConstraintComponentData& data, const SplitSolution& s,
                       SplitKKTResidual& kkt_residual) const override;

  void condenseSlackAndDual(const ContactStatus& contact_status, 
                            ConstraintComponentData& data, 
                            SplitKKTMatrix& kkt_matrix,
                            SplitKKTResidual& kkt_residual) const override;

  void expandSlackAndDual(const ContactStatus& contact_status, 
                          ConstraintComponentData& data, 
                          const SplitDirection& d) const override; 

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

#endif // ROBOTOC_WRENCH_FRICTION_CONE_HPP_ 