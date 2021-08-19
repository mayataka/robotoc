#ifndef IDOCP_LINE_SEARCH_HPP_
#define IDOCP_LINE_SEARCH_HPP_

#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/utils/aligned_vector.hpp"
#include "idocp/cost/cost_function.hpp"
#include "idocp/constraints/constraints.hpp"
#include "idocp/hybrid/contact_sequence.hpp"
#include "idocp/ocp/ocp.hpp"
#include "idocp/ocp/solution.hpp"
#include "idocp/ocp/direction.hpp"
#include "idocp/ocp/kkt_residual.hpp"
#include "idocp/line_search/line_search_filter.hpp"


namespace idocp {

///
/// @class LineSearch 
/// @brief Line search for optimal control problems.
///
class LineSearch {
public:
  ///
  /// @brief Construct a line search.
  /// @param[in] robot Robot model. 
  /// @param[in] N Number of discretization of the horizon. Must be more than 1. 
  /// @param[in] max_num_impulse Maximum number of the impulse on the horizon. 
  /// Must be non-negative. 
  /// @param[in] nthreads Number of the threads in solving the optimal control 
  /// problem. Must be positive. Default is 1.
  /// @param[in] step_size_reduction_rate Reduction rate of the step size. 
  /// Defalt is 0.75.
  /// @param[in] min_step_size Minimum step size. Default is 0.05.
  ///
  LineSearch(const Robot& robot, const int N, const int max_num_impulse=0, 
             const int nthreads=1, const double step_size_reduction_rate=0.75, 
             const double min_step_size=0.05);

  ///
  /// @brief Default constructor. 
  ///
  LineSearch();

  ///
  /// @brief Destructor. 
  ///
  ~LineSearch();

  ///
  /// @brief Default copy constructor. 
  ///
  LineSearch(const LineSearch&) = default;

  ///
  /// @brief Default copy assign operator. 
  ///
  LineSearch& operator=(const LineSearch&) = default;

  ///
  /// @brief Default move constructor. 
  ///
  LineSearch(LineSearch&&) noexcept = default;

  ///
  /// @brief Default move assign operator. 
  ///
  LineSearch& operator=(LineSearch&&) noexcept = default;

  ///
  /// @brief Compute primal step size by fliter line search method. 
  /// @param[in, out] ocp optimal control problem.
  /// @param[in] robots aligned_vector of Robot.
  /// @param[in] contact_sequence Contact sequence. 
  /// @param[in] q Initial configuration.
  /// @param[in] v Initial generalized velocity.
  /// @param[in] s Solution. 
  /// @param[in] d Direction. 
  /// @param[in] max_primal_step_size Maximum primal step size. 
  ///
  double computeStepSize(OCP& ocp, aligned_vector<Robot>& robots,
                         const ContactSequence& contact_sequence, 
                         const Eigen::VectorXd& q, const Eigen::VectorXd& v, 
                         const Solution& s, const Direction& d,
                         const double max_primal_step_size);

  ///
  /// @brief Clear the line search filter. 
  ///
  void clearFilter();

  ///
  /// @brief Checks wheather the line search filter is empty or not. 
  /// @return true if the filter is empty. false if not.
  ///
  bool isFilterEmpty() const;

private:
  LineSearchFilter filter_;
  int max_num_impulse_, nthreads_;
  double step_size_reduction_rate_, min_step_size_;
  Eigen::VectorXd costs_, costs_impulse_, costs_aux_, costs_lift_, violations_, 
                  violations_impulse_, violations_aux_, violations_lift_; 
  Solution s_trial_;
  KKTResidual kkt_residual_;

  void computeCostAndViolation(OCP& ocp, aligned_vector<Robot>& robots,
                               const ContactSequence& contact_sequence, 
                               const Eigen::VectorXd& q, 
                               const Eigen::VectorXd& v, const Solution& s,
                               const double primal_step_size_for_barrier=0);

  void computeSolutionTrial(const OCP& ocp, const aligned_vector<Robot>& robots, 
                            const Solution& s, const Direction& d, 
                            const double step_size);

  static void computeSolutionTrial(const Robot& robot, const SplitSolution& s, 
                                   const SplitDirection& d, 
                                   const double step_size, 
                                   SplitSolution& s_trial) {
    s_trial.setContactStatus(s);
    robot.integrateConfiguration(s.q, d.dq(), step_size, s_trial.q);
    s_trial.v = s.v + step_size * d.dv();
    s_trial.a = s.a + step_size * d.da();
    s_trial.u = s.u + step_size * d.du;
    if (s.hasActiveContacts()) {
      s_trial.f_stack() = s.f_stack() + step_size * d.df();
      s_trial.set_f_vector();
    }
  }

  static void computeSolutionTrial(const Robot& robot, 
                                   const ImpulseSplitSolution& s, 
                                   const ImpulseSplitDirection& d, 
                                   const double step_size, 
                                   ImpulseSplitSolution& s_trial) {
    s_trial.setImpulseStatus(s);
    robot.integrateConfiguration(s.q, d.dq(), step_size, s_trial.q);
    s_trial.v  = s.v  + step_size * d.dv();
    s_trial.dv = s.dv + step_size * d.ddv();
    s_trial.f_stack() = s.f_stack() + step_size * d.df();
    s_trial.set_f_vector();
  }

  void clearCosts() {
    costs_.setZero();
    costs_impulse_.setZero();
    costs_aux_.setZero();
    costs_lift_.setZero();
  }

  void clearViolations() {
    violations_.setZero();
    violations_impulse_.setZero();
    violations_aux_.setZero();
    violations_lift_.setZero();
  }

  double totalCosts() const {
    return (costs_.sum()+costs_impulse_.sum()+costs_aux_.sum()
            +costs_lift_.sum());
  }

  double totalViolations() const {
    return (violations_.sum()+violations_impulse_.sum()+violations_aux_.sum()
            +violations_lift_.sum());
  }

};

} // namespace idocp 

#endif // IDOCP_LINE_SEARCH_HPP_ 