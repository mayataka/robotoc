#ifndef IDOCP_OCP_HPP_
#define IDOCP_OCP_HPP_ 

#include <vector>
#include <memory>

#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/ocp/split_ocp.hpp"
#include "idocp/ocp/terminal_ocp.hpp"
#include "idocp/ocp/line_search_filter.hpp"
#include "idocp/cost/cost_function_interface.hpp"
#include "idocp/constraints/constraints_interface.hpp"


namespace idocp {

class OCP {
public:
  // Constructor. 
  OCP(const Robot& robot, const std::shared_ptr<CostFunctionInterface>& cost,
      const std::shared_ptr<ConstraintsInterface>& constraints, const double T, 
      const int N, const int num_proc=1);

  ~OCP();

  // Prohibit default copy constructor due to unique_ptr.
  OCP(const OCP&) = delete;

  // Prohibit default copy operator due to unique_ptr.
  OCP& operator=(const OCP&) = delete;

  // Use default move constructor.
  OCP(OCP&&) = default;

  // Use default move operator.
  OCP& operator=(OCP&&) = default;

  void solveLQR(const double t, const Eigen::VectorXd& q, 
                const Eigen::VectorXd& v, const bool use_line_search=true);

  void getInitialControlInput(Eigen::VectorXd& u);

  void getStateFeedbackGain(Eigen::MatrixXd& Kq, Eigen::MatrixXd& Kv);

  void setStateTrajectory(const Eigen::VectorXd& q, const Eigen::VectorXd& v);

  void setStateTrajectory(const Eigen::VectorXd& q0, const Eigen::VectorXd& v0,
                          const Eigen::VectorXd& qN, const Eigen::VectorXd& vN);

  double KKTError(const double t, const Eigen::VectorXd& q, 
                  const Eigen::VectorXd& v);

  void printSolution() const;

private:

  bool isCurrentSolutionFeasible();

  void initConstraints();

  void activateAllContacts();

  std::vector<SplitOCP> split_ocps_;
  TerminalOCP terminal_ocp_;
  std::vector<Robot> robots_;
  LineSearchFilter filter_;
  double T_, dtau_, step_size_reduction_rate_, min_step_size_;
  int N_, num_proc_;
  std::vector<Eigen::VectorXd> q_, v_, a_, u_, beta_, f_, mu_, lmd_, gmm_, 
                               dq_, dv_, sq_, sv_;
  std::vector<Eigen::MatrixXd> Pqq_, Pqv_, Pvq_, Pvv_;
  Eigen::VectorXd primal_step_sizes_, dual_step_sizes_, costs_, 
                  constraints_violations_, cost_derivative_dot_direction_;
  std::vector<std::vector<bool>> contact_sequence_;
};

} // namespace idocp 


#endif // IDOCP_OCP_HPP_