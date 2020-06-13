#ifndef IDOCP_SPLIT_OCP_HPP_
#define IDOCP_SPLIT_OCP_HPP_

#include <vector>

#include "Eigen/Core"

#include "robot/robot.hpp"
#include "cost/cost_function_interface.hpp"
#include "constraints/constraints_interface.hpp"
#include "constraints/joint_space_constraints.hpp"


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
 
  // Copy constructor.
  SplitOCP(const SplitOCP& other) = default;

  // Copy operator.
  SplitOCP& operator=(const SplitOCP& other) = default;

  void initConstraints(Robot& robot, const double dtau,
                       const Eigen::VectorXd& q, 
                       const Eigen::VectorXd& v, const Eigen::VectorXd& a);

  void linearizeOCP(Robot& robot, const double t, const double dtau, 
                    const Eigen::VectorXd& lmd, const Eigen::VectorXd& gmm, 
                    const Eigen::VectorXd& q, const Eigen::VectorXd& v, 
                    const Eigen::VectorXd& a, const Eigen::VectorXd& lmd_next, 
                    const Eigen::VectorXd& gmm_next, 
                    const Eigen::VectorXd& q_next,
                    const Eigen::VectorXd& v_next);

  void linearizeOCP(Robot& robot, const double t, const Eigen::VectorXd& lmd, 
                    const Eigen::VectorXd& gmm, const Eigen::VectorXd& q, 
                    const Eigen::VectorXd& v, Eigen::MatrixXd& Qqq, 
                    Eigen::MatrixXd& Qqv, Eigen::MatrixXd& Qvq, 
                    Eigen::MatrixXd& Qvv, Eigen::VectorXd& Qq, 
                    Eigen::VectorXd& Qv);

  void backwardRecursion(const double dtau, const Eigen::MatrixXd& Pqq_next, 
                         const Eigen::MatrixXd& Pqv_next, 
                         const Eigen::MatrixXd& Pvq_next, 
                         const Eigen::MatrixXd& Pvv_next, 
                         const Eigen::VectorXd& sq_next, 
                         const Eigen::VectorXd& sv_next, Eigen::MatrixXd& Pqq, 
                         Eigen::MatrixXd& Pqv, Eigen::MatrixXd& Pvq, 
                         Eigen::MatrixXd& Pvv, Eigen::VectorXd& sq, 
                         Eigen::VectorXd& sv);

  void forwardRecursion(const double dtau, const Eigen::VectorXd& dq,   
                        const Eigen::VectorXd& dv, Eigen::VectorXd& da, 
                        Eigen::VectorXd& dq_next, 
                        Eigen::VectorXd& dv_next) const;

  double stepLength(Robot& robot, const double dtau, const Eigen::VectorXd& dq, 
                    const Eigen::VectorXd& dv, const Eigen::VectorXd& da, 
                    Eigen::VectorXd& q, Eigen::VectorXd& v, Eigen::VectorXd& a);

  void updateOCP(Robot& robot, const double dtau,   const Eigen::VectorXd& dq, 
                 const Eigen::VectorXd& dv, const Eigen::VectorXd& da, 
                 const Eigen::MatrixXd& Pqq, const Eigen::MatrixXd& Pqv, 
                 const Eigen::MatrixXd& Pvq, const Eigen::MatrixXd& Pvv, 
                 const Eigen::VectorXd& sq, const Eigen::VectorXd& sv, 
                 Eigen::VectorXd& q, Eigen::VectorXd& v, Eigen::VectorXd& a, 
                 Eigen::VectorXd& lmd, Eigen::VectorXd& gmm);

  void updateOCP(Robot& robot, const Eigen::VectorXd& dq, 
                 const Eigen::VectorXd& dv, const Eigen::MatrixXd& Pqq, 
                 const Eigen::MatrixXd& Pqv, const Eigen::MatrixXd& Pvq, 
                 const Eigen::MatrixXd& Pvv, const Eigen::VectorXd& sq, 
                 const Eigen::VectorXd& sv, Eigen::VectorXd& q, 
                 Eigen::VectorXd& v, Eigen::VectorXd& lmd, 
                 Eigen::VectorXd& gmm) const;

  double squaredOCPErrorNorm(Robot& robot, const double t, const double dtau, 
                             const Eigen::VectorXd& lmd, 
                             const Eigen::VectorXd& gmm, 
                             const Eigen::VectorXd& q, const Eigen::VectorXd& v, 
                             const Eigen::VectorXd& a, 
                             const Eigen::VectorXd& lmd_next, 
                             const Eigen::VectorXd& gmm_next, 
                             const Eigen::VectorXd& q_next,
                             const Eigen::VectorXd& v_next);

  double squaredOCPErrorNorm(Robot& robot, const double t, 
                             const Eigen::VectorXd& lmd, 
                             const Eigen::VectorXd& gmm, 
                             const Eigen::VectorXd& q, 
                             const Eigen::VectorXd& v);

private:
  CostFunctionInterface *cost_;
  ConstraintsInterface *constraints_;
  JointSpaceConstraints joint_constraints_;
  unsigned int dimq_, dimv_;
  Eigen::VectorXd u_, lu_, lq_, lv_, la_, k_;
  Eigen::VectorXd q_res_, v_res_, a_res_;
  Eigen::MatrixXd luu_, du_dq_, du_dv_, du_da_, Qqq_, Qqv_, Qqa_, Qvq_, Qvv_, 
                  Qva_, Qaa_, Ginv_, Kq_, Kv_;
};

} // namespace idocp


#endif // IDOCP_SPLIT_OCP_HPP_