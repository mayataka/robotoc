#ifndef IDOCP_CONSTRAINT_CONPONENT_BASE_HPP_
#define IDOCP_CONSTRAINT_CONPONENT_BASE_HPP_

#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/constraints/constraint_component_data.hpp"


namespace idocp {

class ConstraintComponentBase {
public:
  ConstraintComponentBase() {}

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
                               ConstraintComponentData& data, const double dtau, 
                               const Eigen::Ref<const Eigen::VectorXd>& a, 
                               const Eigen::Ref<const Eigen::VectorXd>& f, 
                               const Eigen::Ref<const Eigen::VectorXd>& q, 
                               const Eigen::Ref<const Eigen::VectorXd>& v, 
                               const Eigen::Ref<const Eigen::VectorXd>& u) const = 0;

  virtual void augmentDualResidual(const Robot& robot, 
                                   ConstraintComponentData& data, 
                                   const double dtau, 
                                   Eigen::Ref<Eigen::VectorXd> lu) const = 0;

  virtual void augmentDualResidual(const Robot& robot, 
                                   ConstraintComponentData& data, 
                                   const double dtau, 
                                   Eigen::Ref<Eigen::VectorXd> la,
                                   Eigen::Ref<Eigen::VectorXd> lf, 
                                   Eigen::Ref<Eigen::VectorXd> lq,  
                                   Eigen::Ref<Eigen::VectorXd> lv) const = 0;

  virtual void condenseSlackAndDual(const Robot& robot, 
                                    ConstraintComponentData& data, 
                                    const double dtau, 
                                    const Eigen::Ref<const Eigen::VectorXd>& a, 
                                    const Eigen::Ref<const Eigen::VectorXd>& f, 
                                    const Eigen::Ref<const Eigen::VectorXd>& q, 
                                    const Eigen::Ref<const Eigen::VectorXd>& v, 
                                    Eigen::Ref<Eigen::MatrixXd> Caa,
                                    Eigen::Ref<Eigen::MatrixXd> Cff, 
                                    Eigen::Ref<Eigen::MatrixXd> Cqa,  
                                    Eigen::Ref<Eigen::MatrixXd> Cvv, 
                                    Eigen::Ref<Eigen::VectorXd> la,
                                    Eigen::Ref<Eigen::VectorXd> lf, 
                                    Eigen::Ref<Eigen::VectorXd> lq, 
                                    Eigen::Ref<Eigen::VectorXd> lv) const = 0;

  virtual void condenseSlackAndDual(const Robot& robot, 
                                    ConstraintComponentData& data, 
                                    const double dtau, 
                                    const Eigen::Ref<const Eigen::VectorXd>& u, 
                                    Eigen::Ref<Eigen::MatrixXd> Cuu, 
                                    Eigen::Ref<Eigen::VectorXd> Cu) const = 0;

  virtual void computeSlackAndDualDirection(const Robot& robot, 
                                            ConstraintComponentData& data, 
                                            const double dtau,
                                            const Eigen::Ref<const Eigen::VectorXd>& da,
                                            const Eigen::Ref<const Eigen::VectorXd>& df,
                                            const Eigen::Ref<const Eigen::VectorXd>& dq,
                                            const Eigen::Ref<const Eigen::VectorXd>& dv,
                                            const Eigen::Ref<const Eigen::VectorXd>& du) const = 0;

  virtual double maxSlackStepSize(const ConstraintComponentData& data) const = 0;

  virtual double maxDualStepSize(const ConstraintComponentData& data) const = 0;

  virtual void updateSlack(ConstraintComponentData& data, 
                           const double step_size) const = 0;

  virtual void updateDual(ConstraintComponentData& data, 
                          const double step_size) const = 0;

  virtual double costSlackBarrier(const ConstraintComponentData& data) const = 0;

  virtual double costSlackBarrier(const ConstraintComponentData& data, 
                                  const double step_size) const = 0;

  virtual double residualL1Nrom(const Robot& robot, 
                                ConstraintComponentData& data, 
                                const double dtau, 
                                const Eigen::Ref<const Eigen::VectorXd>& a, 
                                const Eigen::Ref<const Eigen::VectorXd>& f, 
                                const Eigen::Ref<const Eigen::VectorXd>& q, 
                                const Eigen::Ref<const Eigen::VectorXd>& v, 
                                const Eigen::Ref<const Eigen::VectorXd>& u) const = 0;

  virtual double residualSquaredNrom(const Robot& robot, 
                                     ConstraintComponentData& data, 
                                     const double dtau,
                                     const Eigen::Ref<const Eigen::VectorXd>& a, 
                                     const Eigen::Ref<const Eigen::VectorXd>& f, 
                                     const Eigen::Ref<const Eigen::VectorXd>& q, 
                                     const Eigen::Ref<const Eigen::VectorXd>& v, 
                                     const Eigen::Ref<const Eigen::VectorXd>& u) const = 0;
  
  virtual int dim() const = 0;

};

} // namespace idocp

#endif // IDOCP_CONSTRAINT_CONPONENT_BASE_HPP_