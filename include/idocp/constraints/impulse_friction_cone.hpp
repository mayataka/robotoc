#ifndef IDOCP_IMPULSE_FRICTION_CONE_HPP_
#define IDOCP_IMPULSE_FRICTION_CONE_HPP_

#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/impulse/impulse_split_solution.hpp"
#include "idocp/impulse/impulse_split_direction.hpp"
#include "idocp/constraints/impulse_constraint_component_base.hpp"
#include "idocp/constraints/constraint_component_data.hpp"
#include "idocp/impulse/impulse_split_kkt_residual.hpp"
#include "idocp/impulse/impulse_split_kkt_matrix.hpp"


namespace idocp {

///
/// @class ImpulseFrictionCone
/// @brief Constraint on the inner-approximated impulse firction cone.
///
class ImpulseFrictionCone final : public ImpulseConstraintComponentBase {
public:
  ///
  /// @brief Constructor. 
  /// @param[in] robot Robot model.
  /// @param[in] mu Friction coefficient. Must be positive.
  /// @param[in] barrier Barrier parameter. Must be positive. Should be small.
  /// Default is 1.0e-04.
  /// @param[in] fraction_to_boundary_rate Must be larger than 0 and smaller 
  /// than 1. Should be between 0.9 and 0.995. Default is 0.995.
  ///
  ImpulseFrictionCone(const Robot& robot, const double mu, 
                      const double barrier=1.0e-04,
                      const double fraction_to_boundary_rate=0.995);

  ///
  /// @brief Default constructor. 
  ///
  ImpulseFrictionCone();

  ///
  /// @brief Destructor. 
  ///
  ~ImpulseFrictionCone();

  ///
  /// @brief Default copy constructor. 
  ///
  ImpulseFrictionCone(const ImpulseFrictionCone&) = default;

  ///
  /// @brief Default copy operator. 
  ///
  ImpulseFrictionCone& operator=(const ImpulseFrictionCone&) = default;

  ///
  /// @brief Default move constructor. 
  ///
  ImpulseFrictionCone(ImpulseFrictionCone&&) noexcept = default;

  ///
  /// @brief Default move assign operator. 
  ///
  ImpulseFrictionCone& operator=(ImpulseFrictionCone&&) noexcept = default;

  ///
  /// @brief Sets the friction coefficient. 
  /// @param[in] mu Friction coefficient. Must be positive.
  ///
  void setFrictionCoefficient(const double mu);

  KinematicsLevel kinematicsLevel() const override;

  void allocateExtraData(ConstraintComponentData& data) const override;

  bool isFeasible(Robot& robot, ConstraintComponentData& data, 
                  const ImpulseSplitSolution& s) const override;

  void setSlackAndDual(Robot& robot, ConstraintComponentData& data, 
                       const ImpulseSplitSolution& s) const override;

  void augmentDualResidual(Robot& robot, ConstraintComponentData& data, 
                           const ImpulseSplitSolution& s,
                           ImpulseSplitKKTResidual& kkt_residual) const override;

  void condenseSlackAndDual(Robot& robot, ConstraintComponentData& data, 
                            const ImpulseSplitSolution& s,
                            ImpulseSplitKKTMatrix& kkt_matrix,
                            ImpulseSplitKKTResidual& kkt_residual) const override;

  void expandSlackAndDual(ConstraintComponentData& data, 
                          const ImpulseSplitSolution& s,
                          const ImpulseSplitDirection& d) const override; 

  void computePrimalAndDualResidual(Robot& robot, ConstraintComponentData& data, 
                                    const ImpulseSplitSolution& s) const override;
  
  int dimc() const override;

  ///
  /// @brief Transforms the contact force from the local coordinate to the 
  /// world coordinate.
  /// @param[in] robot Robot model. Kinematics must be updated.
  /// @param[in] contact_frame_id Index of the contact frame.
  /// @param[in] f_local Contact force expressed in the local frame.
  /// @param[out] f_world Contact force expressed in the world frame. Size must 
  /// be 3.
  ///
  template <typename VectorType>
  static void fLocal2World(const Robot& robot, const int contact_frame_id, 
                           const Eigen::Vector3d& f_local,
                           const Eigen::MatrixBase<VectorType>& f_world) {
    assert(f_world.size() == 3);
    const_cast<Eigen::MatrixBase<VectorType>&>(f_world).noalias()
        = robot.frameRotation(contact_frame_id) * f_local;
  }

  ///
  /// @brief Computes the friction cone residual.
  /// @param[in] mu Friction coefficient. Must be positive.
  /// @param[in] f Contact force expressed in the world frame. Size must be 3.
  /// @param[out] res Friction cone residual. Size must be 5.
  ///
  template <typename VectorType1, typename VectorType2>
  static void frictionConeResidual(const double mu, 
                                   const Eigen::MatrixBase<VectorType1>& f,
                                   const Eigen::MatrixBase<VectorType2>& res) {
    assert(mu > 0);
    assert(f.size() == 3);
    assert(res.size() == 5);
    const_cast<Eigen::MatrixBase<VectorType2>&>(res).coeffRef(0) = - f.coeff(2);
    const_cast<Eigen::MatrixBase<VectorType2>&>(res).coeffRef(1) 
        = f.coeff(0) - mu * f.coeff(2) / std::sqrt(2);
    const_cast<Eigen::MatrixBase<VectorType2>&>(res).coeffRef(2) 
        = - f.coeff(0) - mu * f.coeff(2) / std::sqrt(2);
    const_cast<Eigen::MatrixBase<VectorType2>&>(res).coeffRef(3) 
        = f.coeff(1) - mu * f.coeff(2) / std::sqrt(2);
    const_cast<Eigen::MatrixBase<VectorType2>&>(res).coeffRef(4) 
        = - f.coeff(1) - mu * f.coeff(2) / std::sqrt(2);
  }

  EIGEN_MAKE_ALIGNED_OPERATOR_NEW 

private:
  int dimv_, dimc_, max_point_contacts_;
  std::vector<int> contact_frame_;
  double mu_;
  Eigen::MatrixXd cone_;
  // Eigen::Matrix<double, 5, 3> cone_;

  Eigen::VectorXd& fW(ConstraintComponentData& data, 
                      const int contact_idx) const {
    return data.r[contact_idx];
  }

  Eigen::VectorXd& r(ConstraintComponentData& data, 
                     const int contact_idx) const {
    return data.r[max_point_contacts_+contact_idx];
  }

  Eigen::MatrixXd& dg_dq(ConstraintComponentData& data, 
                         const int contact_idx) const {
    return data.J[contact_idx];
  }

  Eigen::MatrixXd& dg_df(ConstraintComponentData& data, 
                         const int contact_idx) const {
    return data.J[max_point_contacts_+contact_idx];
  }

  Eigen::MatrixXd& dfW_dq(ConstraintComponentData& data, 
                          const int contact_idx) const {
    return data.J[2*max_point_contacts_+contact_idx];
  }

  Eigen::MatrixXd& r_dg_df(ConstraintComponentData& data, 
                           const int contact_idx) const {
    return data.J[3*max_point_contacts_+contact_idx];
  }

};

} // namespace idocp

#endif // IDOCP_IMPULSE_FRICTION_CONE_HPP_ 