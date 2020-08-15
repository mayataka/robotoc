#ifndef IDOCP_CONSTRAINT_COMPONENT_BASE_HPP_
#define IDOCP_CONSTRAINT_COMPONENT_BASE_HPP_

#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/ocp/split_solution.hpp"
#include "idocp/ocp/split_direction.hpp"
#include "idocp/constraints/constraint_component_data.hpp"


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

  virtual bool isFeasible(const Robot& robot, ConstraintComponentData& data, 
                          const Eigen::Ref<const Eigen::VectorXd>& a, 
                          const Eigen::Ref<const Eigen::VectorXd>& f, 
                          const Eigen::Ref<const Eigen::VectorXd>& q, 
                          const Eigen::Ref<const Eigen::VectorXd>& v,
                          const Eigen::Ref<const Eigen::VectorXd>& u) const = 0;

  virtual void setSlackAndDual(const Robot& robot, 
                               ConstraintComponentData& data, 
                               const double dtau, 
                               const Eigen::Ref<const Eigen::VectorXd>& a, 
                               const Eigen::Ref<const Eigen::VectorXd>& f, 
                               const Eigen::Ref<const Eigen::VectorXd>& q, 
                               const Eigen::Ref<const Eigen::VectorXd>& v,
                               const Eigen::Ref<const Eigen::VectorXd>& u) const = 0;

  virtual void augmentDualResidual(const Robot& robot, 
                                   ConstraintComponentData& datas,
                                   const double dtau,
                                   Eigen::Ref<Eigen::VectorXd> la, 
                                   Eigen::Ref<Eigen::VectorXd> lf,
                                   Eigen::Ref<Eigen::VectorXd> lq, 
                                   Eigen::Ref<Eigen::VectorXd> lv) const = 0;

  virtual void augmentDualResidual(const Robot& robot,
                                   ConstraintComponentData& datas,
                                   const double dtau,
                                   Eigen::Ref<Eigen::VectorXd> lu) const = 0;

  virtual void condenseSlackAndDual(const Robot& robot, 
                                    ConstraintComponentData& datas,
                                    const double dtau, 
                                    const Eigen::Ref<const Eigen::VectorXd>& a, 
                                    const Eigen::Ref<const Eigen::VectorXd>& f, 
                                    const Eigen::Ref<const Eigen::VectorXd>& q, 
                                    const Eigen::Ref<const Eigen::VectorXd>& v,
                                    Eigen::Ref<Eigen::MatrixXd> Qaa, 
                                    Eigen::Ref<Eigen::MatrixXd> Qff,
                                    Eigen::Ref<Eigen::MatrixXd> Qqq, 
                                    Eigen::Ref<Eigen::MatrixXd> Qvv,
                                    Eigen::Ref<Eigen::VectorXd> la, 
                                    Eigen::Ref<Eigen::VectorXd> lf,
                                    Eigen::Ref<Eigen::VectorXd> lq, 
                                    Eigen::Ref<Eigen::VectorXd> lv) const = 0;

  virtual void condenseSlackAndDual(const Robot& robot, 
                                    ConstraintComponentData& data, 
                                    const double dtau, 
                                    const Eigen::Ref<const Eigen::VectorXd>& u, 
                                    Eigen::Ref<Eigen::MatrixXd> Quu,
                                    Eigen::Ref<Eigen::VectorXd> lu) const = 0;

  virtual void computeSlackAndDualDirection(
      const Robot& robot, ConstraintComponentData& data, const double dtau, 
      const Eigen::Ref<const Eigen::VectorXd>& da, 
      const Eigen::Ref<const Eigen::VectorXd>& df, 
      const Eigen::Ref<const Eigen::VectorXd>& dq, 
      const Eigen::Ref<const Eigen::VectorXd>& dv, 
      const Eigen::Ref<const Eigen::VectorXd>& du) const = 0;

  virtual double residualL1Nrom(
      const Robot& robot, ConstraintComponentData& data, 
      const double dtau, const Eigen::Ref<const Eigen::VectorXd>& a, 
      const Eigen::Ref<const Eigen::VectorXd>& f, 
      const Eigen::Ref<const Eigen::VectorXd>& q, 
      const Eigen::Ref<const Eigen::VectorXd>& v, 
      const Eigen::Ref<const Eigen::VectorXd>& u) const = 0;

  virtual double residualSquaredNrom(
      const Robot& robot, ConstraintComponentData& data, 
      const double dtau, const Eigen::Ref<const Eigen::VectorXd>& a, 
      const Eigen::Ref<const Eigen::VectorXd>& f, 
      const Eigen::Ref<const Eigen::VectorXd>& q, 
      const Eigen::Ref<const Eigen::VectorXd>& v, 
      const Eigen::Ref<const Eigen::VectorXd>& u) const = 0;
  
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
  virtual void setSlackAndDualPositive(
      Eigen::Ref<Eigen::VectorXd> slack, 
      Eigen::Ref<Eigen::VectorXd> dual) const final;

  virtual void computeDualityResidual(
      const Eigen::Ref<const Eigen::VectorXd>& slack, 
      const Eigen::Ref<const Eigen::VectorXd>& dual, 
      Eigen::Ref<Eigen::VectorXd> duality) const final;

  virtual void computeDualDirection(
      const Eigen::Ref<const Eigen::VectorXd>& slack, 
      const Eigen::Ref<const Eigen::VectorXd>& dslack, 
      const Eigen::Ref<const Eigen::VectorXd>& dual, 
      const Eigen::Ref<const Eigen::VectorXd>& duality, 
      Eigen::Ref<Eigen::VectorXd> ddual) const final;

  virtual double fractionToBoundary(
      const Eigen::Ref<const Eigen::VectorXd>& vec, 
      const Eigen::Ref<const Eigen::VectorXd>& dvec) const final;

private:
  double barrier_, fraction_to_boundary_rate_;

};

} // namespace idocp


#endif // IDOCP_CONSTRAINT_COMPONENT_BASE_HPP_