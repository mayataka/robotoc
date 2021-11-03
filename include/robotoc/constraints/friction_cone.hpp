#ifndef ROBOTOC_FRICTION_CONE_HPP_
#define ROBOTOC_FRICTION_CONE_HPP_

#include "Eigen/Core"

#include "robotoc/robot/robot.hpp"
#include "robotoc/ocp/split_solution.hpp"
#include "robotoc/ocp/split_direction.hpp"
#include "robotoc/constraints/constraint_component_base.hpp"
#include "robotoc/constraints/constraint_component_data.hpp"
#include "robotoc/ocp/split_kkt_residual.hpp"
#include "robotoc/ocp/split_kkt_matrix.hpp"


namespace robotoc {

///
/// @class FrictionCone
/// @brief Constraint on the inner-approximated firction cone.
///
class FrictionCone final : public ConstraintComponentBase {
public:
  using Vector5d = Eigen::Matrix<double, 5, 1>;

  ///
  /// @brief Constructor. 
  /// @param[in] robot Robot model.
  /// @param[in] mu Friction coefficient. Must be positive.
  /// @param[in] barrier Barrier parameter. Must be positive. Should be small.
  /// Default is 1.0e-04.
  /// @param[in] fraction_to_boundary_rule Parameter of the 
  /// fraction-to-boundary-rule Must be larger than 0 and smaller than 1. 
  /// Should be between 0.9 and 0.995. Default is 0.995.
  ///
  FrictionCone(const Robot& robot, const double mu, 
               const double barrier=1.0e-04,
               const double fraction_to_boundary_rule=0.995);

  ///
  /// @brief Default constructor. 
  ///
  FrictionCone();

  ///
  /// @brief Destructor. 
  ///
  ~FrictionCone();

  ///
  /// @brief Default copy constructor. 
  ///
  FrictionCone(const FrictionCone&) = default;

  ///
  /// @brief Default copy operator. 
  ///
  FrictionCone& operator=(const FrictionCone&) = default;

  ///
  /// @brief Default move constructor. 
  ///
  FrictionCone(FrictionCone&&) noexcept = default;

  ///
  /// @brief Default move assign operator. 
  ///
  FrictionCone& operator=(FrictionCone&&) noexcept = default;

  ///
  /// @brief Sets the friction coefficient. 
  /// @param[in] mu Friction coefficient. Must be positive.
  ///
  void setFrictionCoefficient(const double mu);

  bool useKinematics() const override;

  KinematicsLevel kinematicsLevel() const override;

  void allocateExtraData(ConstraintComponentData& data) const override;

  bool isFeasible(Robot& robot, ConstraintComponentData& data, 
                  const SplitSolution& s) const override;

  void setSlack(Robot& robot, ConstraintComponentData& data, 
                const SplitSolution& s) const override;

  void evalConstraint(Robot& robot, ConstraintComponentData& data, 
                      const SplitSolution& s) const override;

  void evalDerivatives(Robot& robot, ConstraintComponentData& data, 
                       const SplitSolution& s,
                       SplitKKTResidual& kkt_residual) const override;

  void condenseSlackAndDual(ConstraintComponentData& data, 
                            const SplitSolution& s, SplitKKTMatrix& kkt_matrix,
                            SplitKKTResidual& kkt_residual) const override;

  void expandSlackAndDual(ConstraintComponentData& data, const SplitSolution& s,
                          const SplitDirection& d) const override; 

  int dimc() const override;

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

  static Eigen::VectorXd& fW(ConstraintComponentData& data, 
                             const int contact_idx) {
    return data.r[contact_idx];
  }

  Eigen::VectorXd& r(ConstraintComponentData& data, 
                     const int contact_idx) const {
    return data.r[max_point_contacts_+contact_idx];
  }

  static Eigen::MatrixXd& dg_dq(ConstraintComponentData& data, 
                                const int contact_idx) {
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

} // namespace robotoc

#endif // ROBOTOC_FRICTION_CONE_HPP_ 