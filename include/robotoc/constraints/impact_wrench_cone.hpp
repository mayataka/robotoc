#ifndef ROBOTOC_IMPACT_WRENCH_CONE_HPP_
#define ROBOTOC_IMPACT_WRENCH_CONE_HPP_

#include "Eigen/Core"

#include "robotoc/robot/robot.hpp"
#include "robotoc/robot/impact_status.hpp"
#include "robotoc/core/split_solution.hpp"
#include "robotoc/core/split_direction.hpp"
#include "robotoc/core/split_kkt_residual.hpp"
#include "robotoc/core/split_kkt_matrix.hpp"
#include "robotoc/constraints/impact_constraint_component_base.hpp"
#include "robotoc/constraints/constraint_component_data.hpp"


namespace robotoc {

///
/// @class ImpactWrenchCone
/// @brief Constraint on the wrench firction cone for surface contacts.
///
class ImpactWrenchCone final : public ImpactConstraintComponentBase {
public:
  ///
  /// @brief Constructor. 
  /// @param[in] robot Robot model.
  /// @param[in] X A length of the rectangular. Must be positive.
  /// @param[in] Y A length of the rectangular. Must be positive.
  ///
  ImpactWrenchCone(const Robot& robot, const double X, const double Y);

  ///
  /// @brief Default constructor. 
  ///
  ImpactWrenchCone();

  ///
  /// @brief Destructor. 
  ///
  ~ImpactWrenchCone();

  ///
  /// @brief Default copy constructor. 
  ///
  ImpactWrenchCone(const ImpactWrenchCone&) = default;

  ///
  /// @brief Default copy operator. 
  ///
  ImpactWrenchCone& operator=(const ImpactWrenchCone&) = default;

  ///
  /// @brief Default move constructor. 
  ///
  ImpactWrenchCone(ImpactWrenchCone&&) noexcept = default;

  ///
  /// @brief Default move assign operator. 
  ///
  ImpactWrenchCone& operator=(ImpactWrenchCone&&) noexcept = default;

  ///
  /// @param[in] X A length of the rectangular. Must be positive.
  /// @param[in] Y A length of the rectangular. Must be positive.
  ///
  void setRectangular(const double X, const double Y);

  KinematicsLevel kinematicsLevel() const override;

  void allocateExtraData(ConstraintComponentData& data) const override;

  bool isFeasible(Robot& robot, const ImpactStatus& impact_status, 
                  ConstraintComponentData& data, 
                  const SplitSolution& s) const override;

  void setSlack(Robot& robot, const ImpactStatus& impact_status, 
                ConstraintComponentData& data, 
                const SplitSolution& s) const override;

  void evalConstraint(Robot& robot, const ImpactStatus& impact_status, 
                      ConstraintComponentData& data, 
                      const SplitSolution& s) const override;

  void evalDerivatives(Robot& robot, const ImpactStatus& impact_status, 
                       ConstraintComponentData& data, 
                       const SplitSolution& s,
                       SplitKKTResidual& kkt_residual) const override;

  void condenseSlackAndDual(const ImpactStatus& impact_status,
                            ConstraintComponentData& data, 
                            SplitKKTMatrix& kkt_matrix,
                            SplitKKTResidual& kkt_residual) const override;

  void expandSlackAndDual(const ImpactStatus& impact_status, 
                          ConstraintComponentData& data, 
                          const SplitDirection& d) const override; 

  int dimc() const override;

  EIGEN_MAKE_ALIGNED_OPERATOR_NEW 

private:
  int dimv_, dimc_, max_num_contacts_;
  std::vector<int> contact_frame_;
  std::vector<ContactType> contact_types_;
  double X_, Y_;

  void computeCone(const double mu, Eigen::MatrixXd& cone) const;
  void updateCone(const double mu, Eigen::MatrixXd& cone) const;

};

} // namespace robotoc

#endif // ROBOTOC_IMPACT_WRENCH_CONE_HPP_