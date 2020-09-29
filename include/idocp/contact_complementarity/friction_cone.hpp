#ifndef IDOCP_FRICTION_CONE_HPP_
#define IDOCP_FRICTION_CONE_HPP_

#include <vector>

#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/ocp/split_solution.hpp"
#include "idocp/ocp/split_direction.hpp"
#include "idocp/contact_complementarity/contact_complementarity_component_base.hpp"
#include "idocp/constraints/constraint_component_data.hpp"
#include "idocp/ocp/kkt_residual.hpp"
#include "idocp/ocp/kkt_matrix.hpp"


namespace idocp {

class FrictionCone final : public ContactComplementarityComponentBase<FrictionCone> {
public:
  FrictionCone(const Robot& robot, const double barrier=1.0e-04, 
                     const double fraction_to_boundary_rate=0.995);

  FrictionCone();

  ~FrictionCone();

  // Use default copy constructor.
  FrictionCone(const FrictionCone&) = default;

  // Use default copy coperator.
  FrictionCone& operator=(const FrictionCone&) = default;

  // Use default move constructor.
  FrictionCone(FrictionCone&&) noexcept = default;

  // Use default move assign coperator.
  FrictionCone& operator=(FrictionCone&&) noexcept = default;

  bool isFeasible_impl(Robot& robot, ConstraintComponentData& data, 
                       const SplitSolution& s) const;

  void setSlackAndDual_impl(Robot& robot, ConstraintComponentData& data, 
                            const double dtau, const SplitSolution& s) const;

  void augmentDualResidual_impl(Robot& robot, ConstraintComponentData& data, 
                                const double dtau, const SplitSolution& s,
                                KKTResidual& kkt_residual);

  void condenseSlackAndDual_impl(Robot& robot, ConstraintComponentData& data, 
                                 const double dtau, const SplitSolution& s,
                                 KKTMatrix& kkt_matrix,
                                 KKTResidual& kkt_residual) const;

  void computeSlackAndDualDirection_impl(Robot& robot, 
                                         ConstraintComponentData& data, 
                                         const double dtau, 
                                         const SplitSolution& s,
                                         const SplitDirection& d) const; 

  double residualL1Nrom_impl(Robot& robot, ConstraintComponentData& data, 
                             const double dtau, const SplitSolution& s) const;

  double squaredKKTErrorNorm_impl(Robot& robot, ConstraintComponentData& data, 
                                  const double dtau, 
                                  const SplitSolution& s) const;

  int dimc_impl() const;

  double maxSlackStepSize_impl(
      const ConstraintComponentData& data, 
      const std::vector<bool>& is_contact_active) const;

  double maxDualStepSize_impl(const ConstraintComponentData& data, 
                              const std::vector<bool>& is_contact_active) const;

  double updateSlack_impl(ConstraintComponentData& data, 
                          const std::vector<bool>& is_contact_active,
                          const double step_size) const;

  double updateDual_impl(ConstraintComponentData& data,
                         const std::vector<bool>& is_contact_active,
                         const double step_size) const;

  double costSlackBarrier_impl(
      const ConstraintComponentData& data, 
      const std::vector<bool>& is_contact_active) const;

  double costSlackBarrier_impl(const ConstraintComponentData& data, 
                               const std::vector<bool>& is_contact_active,
                               const double step_size) const;
 
  void setFrictionCoefficient(const Robot& robot);

  static double frictionConeResidual(const double mu, const Eigen::Vector3d& f);

  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

private:
  int dimc_;
  std::vector<double> mu_;
  std::vector<Eigen::Vector3d> friction_cone_derivative_;

};

} // namespace idocp

#include "idocp/contact_complementarity/friction_cone.hxx"

#endif // IDOCP_FRICTION_CONE_HPP_ 