#ifndef IDOCP_SPLIT_OCP_HPP_
#define IDOCP_SPLIT_OCP_HPP_

#include <vector>
#include <utility>

#include "Eigen/Core"

#include "robot/robot.hpp"
#include "cost/cost_function_interface.hpp"
#include "constraints/constraints_interface.hpp"
#include "constraints/pdipm/joint_space_constraints_pdipm.hpp"
#include "constraints/fb/joint_space_constraints_fb.hpp"


namespace idocp {

class SplitOCP {
public:
  // Constructor. Does not allocate raw arrays.
  // Argments:
  //    model: The pinocchio model. Before call this function, pinocchio model
  //      must be initialized by pinocchio::buildModel().
  SplitOCP(const Robot& robot, const CostFunctionInterface* cost,
           const ConstraintsInterface* constraints);

  ~SplitOCP();
 
  // Use default copy constructor.
  SplitOCP(const SplitOCP& other) = default;

  // Use default copy operator.
  SplitOCP& operator=(const SplitOCP& other) = default;

  bool isFeasible(Robot& robot, const Eigen::VectorXd& q, 
                  const Eigen::VectorXd& v, const Eigen::VectorXd& a, 
                  const Eigen::VectorXd& u);

  void initConstraints(Robot& robot, const double dtau,
                       const Eigen::VectorXd& q, const Eigen::VectorXd& v, 
                       const Eigen::VectorXd& a, const Eigen::VectorXd& u);

  void linearizeOCP(Robot& robot, const double t, const double dtau, 
                    const Eigen::VectorXd& lmd, const Eigen::VectorXd& gmm, 
                    const Eigen::VectorXd& q, const Eigen::VectorXd& v, 
                    const Eigen::VectorXd& a, const Eigen::VectorXd& u, 
                    const Eigen::VectorXd& lmd_next, 
                    const Eigen::VectorXd& gmm_next, 
                    const Eigen::VectorXd& q_next,
                    const Eigen::VectorXd& v_next);

  void linearizeOCP(Robot& robot, const double t, const Eigen::VectorXd& lmd, 
                    const Eigen::VectorXd& gmm, const Eigen::VectorXd& q, 
                    const Eigen::VectorXd& v, Eigen::MatrixXd& Qqq, 
                    Eigen::MatrixXd& Qqv, Eigen::MatrixXd& Qvq, 
                    Eigen::MatrixXd& Qvv, Eigen::VectorXd& Qq, 
                    Eigen::VectorXd& Qv);

  void backwardRiccatiRecursion(const double dtau, 
                                const Eigen::MatrixXd& Pqq_next, 
                                const Eigen::MatrixXd& Pqv_next, 
                                const Eigen::MatrixXd& Pvq_next, 
                                const Eigen::MatrixXd& Pvv_next, 
                                const Eigen::VectorXd& sq_next, 
                                const Eigen::VectorXd& sv_next, 
                                Eigen::MatrixXd& Pqq, Eigen::MatrixXd& Pqv, 
                                Eigen::MatrixXd& Pvq, Eigen::MatrixXd& Pvv, 
                                Eigen::VectorXd& sq, Eigen::VectorXd& sv);

  void forwardRiccatiRecursion(const double dtau, const Eigen::VectorXd& dq,   
                               const Eigen::VectorXd& dv, 
                               Eigen::VectorXd& dq_next, 
                               Eigen::VectorXd& dv_next);

  void computeCondensedDirection(Robot& robot, const double dtau, 
                                 const Eigen::VectorXd& dq, 
                                 const Eigen::VectorXd& dv);
 
  double maxPrimalStepSize();

  double maxDualStepSize();

  double stageCostDerivativeDotDirection(Robot& robot, const double t, 
                                         const double dtau, 
                                         const Eigen::VectorXd& q, 
                                         const Eigen::VectorXd& v, 
                                         const Eigen::VectorXd& a, 
                                         const Eigen::VectorXd& u,
                                         const Eigen::VectorXd& dq,
                                         const Eigen::VectorXd& dv);

  double terminalCostDerivativeDotDirection(Robot& robot, const double t, 
                                            const Eigen::VectorXd& q, 
                                            const Eigen::VectorXd& v, 
                                            const Eigen::VectorXd& dq,
                                            const Eigen::VectorXd& dv);

  std::pair<double, double> stageCostAndConstraintsViolation(
      Robot& robot, const double t, const double dtau, 
      const Eigen::VectorXd& q, const Eigen::VectorXd& v, 
      const Eigen::VectorXd& a, const Eigen::VectorXd& u);

  std::pair<double, double> stageCostAndConstraintsViolation(
      Robot& robot, const double step_size, const double t, const double dtau, 
      const Eigen::VectorXd& q, const Eigen::VectorXd& v, 
      const Eigen::VectorXd& a, const Eigen::VectorXd& u, 
      const Eigen::VectorXd& q_next, const Eigen::VectorXd& v_next, 
      const Eigen::VectorXd& dq, const Eigen::VectorXd& dv, 
      const Eigen::VectorXd& dq_next, const Eigen::VectorXd& dv_next);

  double terminalCost(Robot& robot, const double t, const Eigen::VectorXd& q, 
                      const Eigen::VectorXd& v);

  double terminalCost(Robot& robot, const double step_size, const double t, 
                      const Eigen::VectorXd& q, const Eigen::VectorXd& v, 
                      const Eigen::VectorXd& dq, const Eigen::VectorXd& dv);

  void updateDual(const double step_size);

  void updatePrimal(Robot& robot, const double step_size, const double dtau, 
                    const Eigen::VectorXd& dq, const Eigen::VectorXd& dv, 
                    const Eigen::MatrixXd& Pqq, const Eigen::MatrixXd& Pqv, 
                    const Eigen::MatrixXd& Pvq, const Eigen::MatrixXd& Pvv, 
                    const Eigen::VectorXd& sq, const Eigen::VectorXd& sv, 
                    Eigen::VectorXd& q, Eigen::VectorXd& v, Eigen::VectorXd& a, 
                    Eigen::VectorXd& u, Eigen::VectorXd& beta, 
                    Eigen::VectorXd& lmd, Eigen::VectorXd& gmm);

  void updatePrimal(Robot& robot, const double step_size, 
                    const Eigen::VectorXd& dq, const Eigen::VectorXd& dv, 
                    const Eigen::MatrixXd& Pqq, const Eigen::MatrixXd& Pqv, 
                    const Eigen::MatrixXd& Pvq, const Eigen::MatrixXd& Pvv, 
                    const Eigen::VectorXd& sq, const Eigen::VectorXd& sv, 
                    Eigen::VectorXd& q, Eigen::VectorXd& v, 
                    Eigen::VectorXd& lmd, Eigen::VectorXd& gmm) const;

  void getStateFeedbackGain(Eigen::MatrixXd& Kq, Eigen::MatrixXd& Kv) const;

  double squaredKKTErrorNorm(Robot& robot, const double t, const double dtau, 
                             const Eigen::VectorXd& lmd, 
                             const Eigen::VectorXd& gmm, 
                             const Eigen::VectorXd& q, const Eigen::VectorXd& v, 
                             const Eigen::VectorXd& a, const Eigen::VectorXd& u, 
                             const Eigen::VectorXd& beta, 
                             const Eigen::VectorXd& lmd_next, 
                             const Eigen::VectorXd& gmm_next, 
                             const Eigen::VectorXd& q_next,
                             const Eigen::VectorXd& v_next);

  double squaredKKTErrorNorm(Robot& robot, const double t, 
                             const Eigen::VectorXd& lmd, 
                             const Eigen::VectorXd& gmm, 
                             const Eigen::VectorXd& q, 
                             const Eigen::VectorXd& v);

private:
  CostFunctionInterface *cost_;
  ConstraintsInterface *constraints_;
  pdipm::JointSpaceConstraints joint_constraints_;
  unsigned int dimq_, dimv_;
  Eigen::VectorXd lq_, lv_, la_, lu_, lu_condensed_, ka_, da_;
  Eigen::VectorXd q_res_, v_res_, a_res_, u_res_, du_;
  Eigen::MatrixXd luu_, du_dq_, du_dv_, du_da_, Qqq_, Qqv_, Qqa_, Qvq_, Qvv_, 
                  Qva_, Qaa_, Ginv_, Kq_, Kv_;
  Eigen::VectorXd q_tmp_, v_tmp_, a_tmp_, u_tmp_, u_res_tmp_;
};

} // namespace idocp


#endif // IDOCP_SPLIT_OCP_HPP_