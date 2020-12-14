#ifndef IDOCP_CONTACT_NORMAL_FORCE_HPP_
#define IDOCP_CONTACT_NORMAL_FORCE_HPP_

#include "idocp/robot/robot.hpp"
#include "idocp/ocp/split_solution.hpp"
#include "idocp/ocp/split_direction.hpp"
#include "idocp/constraints/constraint_component_base.hpp"
#include "idocp/constraints/constraint_component_data.hpp"
#include "idocp/ocp/split_kkt_residual.hpp"
#include "idocp/ocp/split_kkt_matrix.hpp"


namespace idocp {

class ContactNormalForce final : public ConstraintComponentBase {
public:
  ContactNormalForce(const Robot& robot, const double barrier=1.0e-04,
                     const double fraction_to_boundary_rate=0.995);

  ContactNormalForce();

  ~ContactNormalForce();

  // Use default copy constructor.
  ContactNormalForce(const ContactNormalForce&) = default;

  // Use default copy coperator.
  ContactNormalForce& operator=(const ContactNormalForce&) = default;

  // Use default move constructor.
  ContactNormalForce(ContactNormalForce&&) noexcept = default;

  // Use default move assign coperator.
  ContactNormalForce& operator=(ContactNormalForce&&) noexcept = default;

  bool useKinematics() const override;

  KinematicsLevel kinematicsLevel() const override;

  bool isFeasible(Robot& robot, ConstraintComponentData& data, 
                  const SplitSolution& s) const override;

  void setSlackAndDual(Robot& robot, ConstraintComponentData& data, 
                       const double dtau, const SplitSolution& s) const override;

  void augmentDualResidual(Robot& robot, ConstraintComponentData& data, 
                           const double dtau, const SplitSolution& s,
                           SplitKKTResidual& kkt_residual) const override;

  void condenseSlackAndDual(Robot& robot, ConstraintComponentData& data, 
                            const double dtau, const SplitSolution& s,
                            SplitKKTMatrix& kkt_matrix,
                            SplitKKTResidual& kkt_residual) const override;

  void computeSlackAndDualDirection(Robot& robot, ConstraintComponentData& data, 
                                    const double dtau, const SplitSolution& s,
                                    const SplitDirection& d) const override; 

  void computePrimalAndDualResidual(Robot& robot, ConstraintComponentData& data, 
                                    const double dtau, 
                                    const SplitSolution& s) const override;
  
  int dimc() const override;

private:
  int dimc_;

};

} // namespace idocp

#endif // IDOCP_CONTACT_NORMAL_FORCE_HPP_ 