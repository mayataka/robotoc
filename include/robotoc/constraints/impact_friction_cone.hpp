#ifndef ROBOTOC_IMPACT_FRICTION_CONE_HPP_
#define ROBOTOC_IMPACT_FRICTION_CONE_HPP_

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
/// @class ImpactFrictionCone
/// @brief Constraint on the inner-approximated impact firction cone.
///
class ImpactFrictionCone final : public ImpactConstraintComponentBase {
public:
  using Vector5d = Eigen::Matrix<double, 5, 1>;

  ///
  /// @brief Constructor. 
  /// @param[in] robot Robot model.
  ///
  ImpactFrictionCone(const Robot& robot);

  ///
  /// @brief Default constructor. 
  ///
  ImpactFrictionCone();

  ///
  /// @brief Destructor. 
  ///
  ~ImpactFrictionCone();

  ///
  /// @brief Default copy constructor. 
  ///
  ImpactFrictionCone(const ImpactFrictionCone&) = default;

  ///
  /// @brief Default copy operator. 
  ///
  ImpactFrictionCone& operator=(const ImpactFrictionCone&) = default;

  ///
  /// @brief Default move constructor. 
  ///
  ImpactFrictionCone(ImpactFrictionCone&&) noexcept = default;

  ///
  /// @brief Default move assign operator. 
  ///
  ImpactFrictionCone& operator=(ImpactFrictionCone&&) noexcept = default;

  KinematicsLevel kinematicsLevel() const override;

  void allocateExtraData(ConstraintComponentData& data) const override;

  bool isFeasible(Robot& robot, const ImpactStatus& impact_status, 
                  const GridInfo& grid_info, const SplitSolution& s,
                  ConstraintComponentData& data) const override;

  void setSlack(Robot& robot, const ImpactStatus& impact_status, 
                const GridInfo& grid_info, const SplitSolution& s,
                ConstraintComponentData& data) const override;

  void evalConstraint(Robot& robot, const ImpactStatus& impact_status, 
                      const GridInfo& grid_info, const SplitSolution& s,
                      ConstraintComponentData& data) const override;

  void evalDerivatives(Robot& robot, const ImpactStatus& impact_status, 
                       const GridInfo& grid_info, const SplitSolution& s,
                       ConstraintComponentData& data,
                       SplitKKTResidual& kkt_residual) const override;

  void condenseSlackAndDual(const ImpactStatus& impact_status, 
                            const GridInfo& grid_info,
                            ConstraintComponentData& data, 
                            SplitKKTMatrix& kkt_matrix,
                            SplitKKTResidual& kkt_residual) const override;

  void expandSlackAndDual(const ImpactStatus& impact_status, 
                          const GridInfo& grid_info,
                          const SplitDirection& d,
                          ConstraintComponentData& data) const override; 

  int dimc() const override;

  ///
  /// @brief Computes the friction cone residual.
  /// @param[in] mu Friction coefficient. Must be positive.
  /// @param[in] f_world Contact force expressed in the world frame. 
  /// Size must be 3.
  /// @param[in] R_surface Rotation matrix of the contact surface.
  /// @param[out] res Friction cone residual. Size must be 5.
  ///
  template <typename VectorType1, typename VectorType2>
  static void frictionConeResidual(const double mu, 
                                   const Eigen::MatrixBase<VectorType1>& f_world,
                                   const Eigen::Matrix3d& R_surface,
                                   const Eigen::MatrixBase<VectorType2>& res) {
    assert(mu > 0);
    assert(f_world.size() == 3);
    assert((R_surface*R_surface.transpose()).isIdentity());
    assert(res.size() == 5);
    const Eigen::Vector3d f_local = R_surface.transpose() * f_world;
    const_cast<Eigen::MatrixBase<VectorType2>&>(res).coeffRef(0) = - f_local.coeff(2);
    const_cast<Eigen::MatrixBase<VectorType2>&>(res).coeffRef(1) 
        = f_local.coeff(0) - mu * f_local.coeff(2) / std::sqrt(2);
    const_cast<Eigen::MatrixBase<VectorType2>&>(res).coeffRef(2) 
        = - f_local.coeff(0) - mu * f_local.coeff(2) / std::sqrt(2);
    const_cast<Eigen::MatrixBase<VectorType2>&>(res).coeffRef(3) 
        = f_local.coeff(1) - mu * f_local.coeff(2) / std::sqrt(2);
    const_cast<Eigen::MatrixBase<VectorType2>&>(res).coeffRef(4) 
        = - f_local.coeff(1) - mu * f_local.coeff(2) / std::sqrt(2);
  }

  EIGEN_MAKE_ALIGNED_OPERATOR_NEW 

private:
  int dimv_, dimc_, max_num_contacts_;
  std::vector<int> contact_frame_;
  std::vector<ContactType> contact_types_;

  Eigen::VectorXd& fW(ConstraintComponentData& data, 
                      const int contact_idx) const {
    return data.r[contact_idx];
  }

  Eigen::VectorXd& r(ConstraintComponentData& data, 
                     const int contact_idx) const {
    return data.r[max_num_contacts_+contact_idx];
  }

  Eigen::MatrixXd& dg_dq(ConstraintComponentData& data, 
                         const int contact_idx) const {
    return data.J[contact_idx];
  }

  Eigen::MatrixXd& dg_df(ConstraintComponentData& data, 
                         const int contact_idx) const {
    return data.J[max_num_contacts_+contact_idx];
  }

  Eigen::MatrixXd& dfW_dq(ConstraintComponentData& data, 
                          const int contact_idx) const {
    return data.J[2*max_num_contacts_+contact_idx];
  }

  Eigen::MatrixXd& r_dg_df(ConstraintComponentData& data, 
                           const int contact_idx) const {
    return data.J[3*max_num_contacts_+contact_idx];
  }

  Eigen::MatrixXd& cone_local(ConstraintComponentData& data, 
                              const int contact_idx) const {
    return data.J[4*max_num_contacts_+contact_idx];
  }

  Eigen::MatrixXd& cone_world(ConstraintComponentData& data, 
                              const int contact_idx) const {
    return data.J[5*max_num_contacts_+contact_idx];
  }

};

} // namespace robotoc

#endif // ROBOTOC_IMPACT_FRICTION_CONE_HPP_ 