#ifndef IDOCP_CONSTRAINT_COMPONENT_BASE_HPP_
#define IDOCP_CONSTRAINT_COMPONENT_BASE_HPP_

#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/ocp/split_solution.hpp"
#include "idocp/ocp/split_direction.hpp"
#include "idocp/constraints/constraint_component_data.hpp"
#include "idocp/ocp/kkt_residual.hpp"
#include "idocp/ocp/kkt_matrix.hpp"


namespace idocp {

class ConstraintComponentBase {
public:
  ConstraintComponentBase(const double barrier, 
                          const double fraction_to_boundary_rate);

  ConstraintComponentBase();

  virtual ~ConstraintComponentBase() {}

  // Use default copy constructor.
  ConstraintComponentBase(const ConstraintComponentBase&) = default;

  // Use default copy coperator.
  ConstraintComponentBase& operator=(const ConstraintComponentBase&) = default;

  // Use default move constructor.
  ConstraintComponentBase(ConstraintComponentBase&&) noexcept = default;

  // Use default move assign coperator.
  ConstraintComponentBase& operator=(ConstraintComponentBase&&) noexcept 
      = default;

  virtual bool useKinematics() const = 0;

  virtual bool isFeasible(const Robot& robot, ConstraintComponentData& data, 
                          const SplitSolution& s) const = 0;

  virtual void setSlackAndDual(const Robot& robot, 
                               ConstraintComponentData& data, const double dtau, 
                               const SplitSolution& s) const = 0;

  virtual void augmentDualResidual(const Robot& robot, 
                                   ConstraintComponentData& data,
                                   const double dtau, 
                                   KKTResidual& kkt_residual) const = 0;

  virtual void augmentDualResidual(const Robot& robot, 
                                   ConstraintComponentData& data,
                                   const double dtau, 
                                   Eigen::VectorXd& lu) const = 0;

  virtual void condenseSlackAndDual(const Robot& robot, 
                                    ConstraintComponentData& data,
                                    const double dtau, const SplitSolution& s, 
                                    KKTMatrix& kkt_matrix,
                                    KKTResidual& kkt_residual) const = 0;

  virtual void condenseSlackAndDual(const Robot& robot, 
                                    ConstraintComponentData& data,
                                    const double dtau, const Eigen::VectorXd& u,
                                    Eigen::MatrixXd& Quu, 
                                    Eigen::VectorXd& lu) const = 0;

  virtual void computeSlackAndDualDirection(
      const Robot& robot, ConstraintComponentData& data, const double dtau, 
      const SplitDirection& d) const = 0;

  virtual double residualL1Nrom(
      const Robot& robot, ConstraintComponentData& data, 
      const double dtau, const SplitSolution& s) const = 0;

  virtual double squaredKKTErrorNorm(
      const Robot& robot, ConstraintComponentData& data, 
      const double dtau, const SplitSolution& s) const = 0;
  
  virtual int dimc() const = 0;

  // Following functions should not be overrided.
  virtual double maxSlackStepSize(
      const ConstraintComponentData& data) const final;

  virtual double maxDualStepSize(
      const ConstraintComponentData& data) const final;

  virtual void updateSlack(ConstraintComponentData& data, 
                           const double step_size) const final;

  virtual void updateDual(ConstraintComponentData& data, 
                          const double step_size) const final;

  virtual double costSlackBarrier(
      const ConstraintComponentData& data) const final;

  virtual double costSlackBarrier(const ConstraintComponentData& data, 
                                  const double step_size) const final;

  virtual void setBarrier(const double barrier) final;

  virtual void setFractionToBoundaryRate(
      const double fraction_to_boundary_rate) final;

protected:
  virtual void setSlackAndDualPositive(Eigen::VectorXd& slack, 
                                       Eigen::VectorXd& dual) const final;

  virtual void computeDuality(const Eigen::VectorXd& slack, 
                              const Eigen::VectorXd& dual, 
                              Eigen::VectorXd& duality) const final;

  virtual void computeDualDirection(const Eigen::VectorXd& slack, 
                                    const Eigen::VectorXd& dual, 
                                    const Eigen::VectorXd& dslack, 
                                    const Eigen::VectorXd& duality, 
                                    Eigen::VectorXd& ddual) const final;

  virtual double fractionToBoundary(const Eigen::VectorXd& vec, 
                                    const Eigen::VectorXd& dvec) const final;

private:
  double barrier_, fraction_to_boundary_rate_;

};

} // namespace idocp

#include "idocp/constraints/constraint_component_base.hxx"

#endif // IDOCP_CONSTRAINT_COMPONENT_BASE_HPP_