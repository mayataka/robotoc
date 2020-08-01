#ifndef IDOCP_QUADRUPED_CONSTRAINTS_HPP_
#define IDOCP_QUADRUPED_CONSTRAINTS_HPP_

#include "Eigen/Core"

#include "robot/robot.hpp"
#include "constraints/constraints_interface.hpp"


namespace idocp {
namespace quadruped {

class Constraints final : public ConstraintsInterface {
public:
  Constraints(const Robot& robot);

  bool isFeasible(const Robot& robot, const Eigen::VectorXd& q, 
                  const Eigen::VectorXd& v, const Eigen::VectorXd& a, 
                  const Eigen::VectorXd& u) override;

  void setSlackAndDual(const Robot& robot, const double dtau, 
                       const Eigen::VectorXd& q, const Eigen::VectorXd& v, 
                       const Eigen::VectorXd& a, 
                       const Eigen::VectorXd& u) override;

  void augmentDualResidual(const Robot& robot, const double dtau, 
                           Eigen::VectorXd& Cu) override;

  void augmentDualResidual(const Robot& robot, const double dtau, 
                           Eigen::VectorXd& Cq, Eigen::VectorXd& Cv, 
                           Eigen::VectorXd& Ca) override;

  void condenseSlackAndDual(const Robot& robot, const double dtau, 
                            const Eigen::VectorXd& q, 
                            const Eigen::VectorXd& v, 
                            const Eigen::VectorXd& a, Eigen::MatrixXd& Cqq, 
                            Eigen::MatrixXd& Cqv, Eigen::MatrixXd& Cqa,  
                            Eigen::MatrixXd& Cvv, Eigen::MatrixXd& Cva,  
                            Eigen::MatrixXd& Caa,  Eigen::VectorXd& Cq, 
                            Eigen::VectorXd& Cv, Eigen::VectorXd& Ca) override;

  void condenseSlackAndDual(const Robot& robot, const double dtau, 
                            const Eigen::VectorXd& u, Eigen::MatrixXd& Cuu, 
                            Eigen::VectorXd& Cu) override;

  void computeSlackAndDualDirection(const Robot& robot, const double dtau,
                                    const Eigen::VectorXd& dq,
                                    const Eigen::VectorXd& dv,
                                    const Eigen::VectorXd& da,
                                    const Eigen::VectorXd& du) override;

  double maxSlackStepSize() override;

  double maxDualStepSize() override;

  void updateSlack(const double step_size) override;

  void updateDual(const double step_size) override;

  double costSlackBarrier() override;

  double costSlackBarrier(const double step_size) override;

  double residualL1Nrom(const Robot& robot, const double dtau,
                        const Eigen::VectorXd& q, const Eigen::VectorXd& v, 
                        const Eigen::VectorXd& a, 
                        const Eigen::VectorXd& u) override;

  double residualSquaredNrom(const Robot& robot, const double dtau,
                             const Eigen::VectorXd& q, const Eigen::VectorXd& v, 
                             const Eigen::VectorXd& a, 
                             const Eigen::VectorXd& u) override;

private:

};

} // namespace quadruped
} // namespace idocp


#endif // IDOCP_QUADRUPED_CONSTRAINTS_HPP_