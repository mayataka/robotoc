#ifndef ROBOTOC_IMPULSE_WRENCH_CONE_HPP_
#define ROBOTOC_IMPULSE_WRENCH_CONE_HPP_

#include "Eigen/Core"

#include "robotoc/robot/robot.hpp"
#include "robotoc/robot/impulse_status.hpp"
#include "robotoc/core/split_solution.hpp"
#include "robotoc/core/split_direction.hpp"
#include "robotoc/core/split_kkt_residual.hpp"
#include "robotoc/core/split_kkt_matrix.hpp"
#include "robotoc/constraints/impulse_constraint_component_base.hpp"
#include "robotoc/constraints/constraint_component_data.hpp"


namespace robotoc {

///
/// @class ImpulseWrenchCone
/// @brief Constraint on the wrench firction cone for surface contacts.
///
class ImpulseWrenchCone final : public ImpulseConstraintComponentBase {
public:
  ///
  /// @brief Constructor. 
  /// @param[in] robot Robot model.
  /// @param[in] X A length of the rectangular. Must be positive.
  /// @param[in] Y A length of the rectangular. Must be positive.
  ///
  ImpulseWrenchCone(const Robot& robot, const double X, const double Y);

  ///
  /// @brief Default constructor. 
  ///
  ImpulseWrenchCone();

  ///
  /// @brief Destructor. 
  ///
  ~ImpulseWrenchCone();

  ///
  /// @brief Default copy constructor. 
  ///
  ImpulseWrenchCone(const ImpulseWrenchCone&) = default;

  ///
  /// @brief Default copy operator. 
  ///
  ImpulseWrenchCone& operator=(const ImpulseWrenchCone&) = default;

  ///
  /// @brief Default move constructor. 
  ///
  ImpulseWrenchCone(ImpulseWrenchCone&&) noexcept = default;

  ///
  /// @brief Default move assign operator. 
  ///
  ImpulseWrenchCone& operator=(ImpulseWrenchCone&&) noexcept = default;

  ///
  /// @param[in] X A length of the rectangular. Must be positive.
  /// @param[in] Y A length of the rectangular. Must be positive.
  ///
  void setRectangular(const double X, const double Y);

  KinematicsLevel kinematicsLevel() const override;

  void allocateExtraData(ConstraintComponentData& data) const override;

  bool isFeasible(Robot& robot, const ImpulseStatus& impulse_status, 
                  ConstraintComponentData& data, 
                  const SplitSolution& s) const override;

  void setSlack(Robot& robot, const ImpulseStatus& impulse_status, 
                ConstraintComponentData& data, 
                const SplitSolution& s) const override;

  void evalConstraint(Robot& robot, const ImpulseStatus& impulse_status, 
                      ConstraintComponentData& data, 
                      const SplitSolution& s) const override;

  void evalDerivatives(Robot& robot, const ImpulseStatus& impulse_status, 
                       ConstraintComponentData& data, 
                       const SplitSolution& s,
                       SplitKKTResidual& kkt_residual) const override;

  void condenseSlackAndDual(const ImpulseStatus& impulse_status,
                            ConstraintComponentData& data, 
                            SplitKKTMatrix& kkt_matrix,
                            SplitKKTResidual& kkt_residual) const override;

  void expandSlackAndDual(const ImpulseStatus& impulse_status, 
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

#endif // ROBOTOC_IMPULSE_WRENCH_CONE_HPP_