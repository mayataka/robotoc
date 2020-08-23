#ifndef IDOCP_OCP_HPP_
#define IDOCP_OCP_HPP_ 

#include <vector>
#include <memory>

#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/ocp/split_ocp.hpp"
#include "idocp/ocp/terminal_ocp.hpp"
#include "idocp/ocp/line_search_filter.hpp"
#include "idocp/ocp/riccati_factorization.hpp"
#include "idocp/cost/cost_function.hpp"
#include "idocp/constraints/constraints.hpp"


namespace idocp {

class OCP {
public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

  // Constructor. 
  OCP(const Robot& robot, const std::shared_ptr<CostFunction>& cost,
      const std::shared_ptr<Constraints>& constraints, const double T, 
      const int N, const int num_proc=1);

  OCP();

  ~OCP();

  // Use default copy constructor.
  OCP(const OCP&) = default;

  // Use default copy assign operator.
  OCP& operator=(const OCP&) = default;

  // Use default move constructor.
  OCP(OCP&&) noexcept = default;

  // Use default move operator.
  OCP& operator=(OCP&&) noexcept = default;

  void updateSolution(const double t, const Eigen::VectorXd& q, 
                      const Eigen::VectorXd& v, const bool use_line_search=true);

  void getInitialControlInput(Eigen::VectorXd& u);

  void getStateFeedbackGain(Eigen::MatrixXd& Kq, Eigen::MatrixXd& Kv);

  bool setStateTrajectory(const Eigen::VectorXd& q, const Eigen::VectorXd& v);

  bool setStateTrajectory(const Eigen::VectorXd& q0, const Eigen::VectorXd& v0,
                          const Eigen::VectorXd& qN, const Eigen::VectorXd& vN);
                          
  void setContactSequence(
      const std::vector<std::vector<bool>>& contact_sequence);

  void resetContactPoint();
  
  void resetLineSearchFilter();

  double KKTError(const double t, const Eigen::VectorXd& q, 
                  const Eigen::VectorXd& v);

  void printSolution() const;

private:

  bool isCurrentSolutionFeasible();

  void initConstraints();

  std::vector<SplitOCP> split_ocps_;
  TerminalOCP terminal_ocp_;
  std::vector<Robot> robots_;
  LineSearchFilter filter_;
  double T_, dtau_, step_size_reduction_rate_, min_step_size_;
  int N_, num_proc_;
  std::vector<SplitSolution> s_;
  std::vector<SplitDirection> d_;
  std::vector<RiccatiFactorization> riccati_;
  Eigen::VectorXd primal_step_sizes_, dual_step_sizes_, costs_, violations_;
  std::vector<std::vector<bool>> contact_sequence_;
};

} // namespace idocp 


#endif // IDOCP_OCP_HPP_