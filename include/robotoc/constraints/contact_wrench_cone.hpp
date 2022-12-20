#ifndef ROBOTOC_CONTACT_WRENCH_CONE_HPP_
#define ROBOTOC_CONTACT_WRENCH_CONE_HPP_

#include "Eigen/Core"

#include "robotoc/robot/robot.hpp"
#include "robotoc/robot/contact_status.hpp"
#include "robotoc/core/split_solution.hpp"
#include "robotoc/core/split_direction.hpp"
#include "robotoc/core/split_kkt_residual.hpp"
#include "robotoc/core/split_kkt_matrix.hpp"
#include "robotoc/constraints/constraint_component_base.hpp"
#include "robotoc/constraints/constraint_component_data.hpp"


namespace robotoc {

///
/// @class ContactWrenchCone
/// @brief Constraint on the contact wrench cone for surface contacts.
///
class ContactWrenchCone final : public ConstraintComponentBase {
public:
  ///
  /// @brief Constructor. 
  /// @param[in] robot Robot model.
  /// @param[in] X A length of the rectangular. Must be positive.
  /// @param[in] Y A length of the rectangular. Must be positive.
  ///
  ContactWrenchCone(const Robot& robot, const double X, const double Y);

  ///
  /// @brief Default constructor. 
  ///
  ContactWrenchCone();

  ///
  /// @brief Destructor. 
  ///
  ~ContactWrenchCone();

  ///
  /// @brief Default copy constructor. 
  ///
  ContactWrenchCone(const ContactWrenchCone&) = default;

  ///
  /// @brief Default copy operator. 
  ///
  ContactWrenchCone& operator=(const ContactWrenchCone&) = default;

  ///
  /// @brief Default move constructor. 
  ///
  ContactWrenchCone(ContactWrenchCone&&) noexcept = default;

  ///
  /// @brief Default move assign operator. 
  ///
  ContactWrenchCone& operator=(ContactWrenchCone&&) noexcept = default;

  ///
  /// @param[in] X A length of the rectangular. Must be positive.
  /// @param[in] Y A length of the rectangular. Must be positive.
  ///
  void setRectangular(const double X, const double Y);

  KinematicsLevel kinematicsLevel() const override;

  void allocateExtraData(ConstraintComponentData& data) const override;

  bool isFeasible(Robot& robot, const ContactStatus& contact_status, 
                  const GridInfo& grid_info, const SplitSolution& s,
                  ConstraintComponentData& data) const override;

  void setSlack(Robot& robot, const ContactStatus& contact_status, 
                const GridInfo& grid_info, const SplitSolution& s,
                ConstraintComponentData& data) const override;

  void evalConstraint(Robot& robot, const ContactStatus& contact_status, 
                      const GridInfo& grid_info, const SplitSolution& s,
                      ConstraintComponentData& data) const override;

  void evalDerivatives(Robot& robot, const ContactStatus& contact_status, 
                       const GridInfo& grid_info, const SplitSolution& s,
                       ConstraintComponentData& data,
                       SplitKKTResidual& kkt_residual) const override;

  void condenseSlackAndDual(const ContactStatus& contact_status, 
                            const GridInfo& grid_info,
                            ConstraintComponentData& data, 
                            SplitKKTMatrix& kkt_matrix,
                            SplitKKTResidual& kkt_residual) const override;

  void expandSlackAndDual(const ContactStatus& contact_status, 
                          const GridInfo& grid_info, const SplitDirection& d, 
                          ConstraintComponentData& data) const override; 

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

#endif // ROBOTOC_CONTACT_WRENCH_CONE_HPP_