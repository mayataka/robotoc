#ifndef IDOCP_PARNMPC_HPP_
#define IDOCP_PARNMPC_HPP_

#include <vector>
#include <memory>

#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/ocp/split_parnmpc.hpp"
#include "idocp/ocp/split_solution.hpp"
#include "idocp/ocp/split_direction.hpp"
#include "idocp/ocp/line_search_filter.hpp"
#include "idocp/cost/cost_function.hpp"
#include "idocp/constraints/constraints.hpp"


namespace idocp {

class ParNMPC {
public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

  // Constructor. 
  ParNMPC(const Robot& robot, const std::shared_ptr<CostFunction>& cost,
          const std::shared_ptr<Constraints>& constraints, 
          const double T, const int N, const int num_proc=1);

  ParNMPC();

  ~ParNMPC();

  // Use default copy constructor.
  ParNMPC(const ParNMPC&) = default;

  // Use default copy assign operator.
  ParNMPC& operator=(const ParNMPC&) = default;

  // Use default move constructor.
  ParNMPC(ParNMPC&&) noexcept = default;

  // Use default move operator.
  ParNMPC& operator=(ParNMPC&&) noexcept = default;

  void updateSolution(const double t, const Eigen::VectorXd& q, 
                      const Eigen::VectorXd& v, 
                      const bool use_line_search=true);

  void getControlInput(const int stage, Eigen::VectorXd& u) const;

  void getStateFeedbackGain(const int stage, Eigen::MatrixXd& Kq, 
                            Eigen::MatrixXd& Kv) const;

  bool setStateTrajectory(const Eigen::VectorXd& q, const Eigen::VectorXd& v);

  bool setStateTrajectory(const Eigen::VectorXd& q0, const Eigen::VectorXd& v0,
                          const Eigen::VectorXd& qN, const Eigen::VectorXd& vN);

  void setAuxiliaryMatrixGuessByTerminalCost(const double t);
                          
  void setContactSequence(
      const std::vector<std::vector<bool>>& contact_sequence);

  void setContactPoint(const std::vector<Eigen::Vector3d>& contact_points);
  
  void clearLineSearchFilter();

  double KKTError(const double t, const Eigen::VectorXd& q, 
                  const Eigen::VectorXd& v);

  void printSolution() const;

private:

  bool isCurrentSolutionFeasible();

  void initConstraints();

  std::vector<SplitParNMPC> split_ocps_;
  std::vector<Robot> robots_;
  LineSearchFilter filter_;
  double T_, dtau_, step_size_reduction_rate_, min_step_size_;
  int N_, num_proc_;
  std::vector<SplitSolution> s_, s_new_;
  std::vector<SplitDirection> d_;
  std::vector<Eigen::MatrixXd> aux_mat_old_;
  Eigen::VectorXd primal_step_sizes_, dual_step_sizes_, costs_, violations_;
  std::vector<std::vector<bool>> contact_sequence_;
};

} // namespace idocp 


#endif // IDOCP_PARNMPC_HPP_