#ifndef IDOCP_LINE_SEARCH_HPP_
#define IDOCP_LINE_SEARCH_HPP_

#include <vector>

#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/utils/aligned_vector.hpp"
#include "idocp/cost/cost_function.hpp"
#include "idocp/constraints/constraints.hpp"
#include "idocp/hybrid/contact_sequence.hpp"
#include "idocp/ocp/ocp.hpp"
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
  /// @param[in] robots std::vector of Robot.
  /// @param[in] contact_sequence Contact sequence. 
  /// @param[in] q Initial configuration.
  /// @param[in] v Initial generalized velocity.
  /// @param[in] s Solution. 
  /// @param[in] d Direction. 
  /// @param[in] max_primal_step_size Maximum primal step size. 
  ///
  template <typename OCPType>
  double computeStepSize(OCPType& ocp, aligned_vector<Robot>& robots,
                         const ContactSequence& contact_sequence, 
                         const Eigen::VectorXd& q, const Eigen::VectorXd& v, 
                         const Solution& s, const Direction& d,
                         const double max_primal_step_size) {
    assert(max_primal_step_size > 0);
    assert(max_primal_step_size <= 1);
    // If filter is empty, augment the current solution to the filter.
    if (filter_.isEmpty()) {
      computeCostAndViolation(ocp, robots, contact_sequence, q, v, s);
      filter_.augment(totalCosts(), totalViolations());
    }
    double primal_step_size = max_primal_step_size;
    while (primal_step_size > min_step_size_) {
      computeSolution(ocp, robots, s, d, primal_step_size);
      computeCostAndViolation(ocp, robots, contact_sequence, q, v, s_try_,
                              primal_step_size);
      const double total_costs = totalCosts();
      const double total_violations = totalViolations();
      if (filter_.isAccepted(total_costs, total_violations)) {
        filter_.augment(total_costs, total_violations);
        break;
      }
      primal_step_size *= step_size_reduction_rate_;
    }
    if (primal_step_size > min_step_size_) {
      return primal_step_size;
    }
    else {
      return min_step_size_;
    }
  }

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
  Solution s_try_;
  KKTResidual kkt_residual_;

  void computeCostAndViolation(OCP& ocp, aligned_vector<Robot>& robots,
                               const ContactSequence& contact_sequence, 
                               const Eigen::VectorXd& q, 
                               const Eigen::VectorXd& v, const Solution& s,
                               const double primal_step_size_for_barrier=0);

  void computeSolution(const OCP& ocp, const aligned_vector<Robot>& robots, 
                          const Solution& s, const Direction& d, 
                          const double step_size);

  static void computeSolution(const Robot& robot, const SplitSolution& s, 
                              const SplitDirection& d, const double step_size, 
                              SplitSolution& s_try) {
    s_try.setContactStatus(s);
    robot.integrateConfiguration(s.q, d.dq(), step_size, s_try.q);
    s_try.v = s.v + step_size * d.dv();
    s_try.a = s.a + step_size * d.da();
    s_try.u = s.u + step_size * d.du;
    if (s.hasActiveContacts()) {
      s_try.f_stack() = s.f_stack() + step_size * d.df();
      s_try.set_f_vector();
    }
  }

  static void computeSolution(const Robot& robot, const ImpulseSplitSolution& s, 
                              const ImpulseSplitDirection& d, 
                              const double step_size, 
                              ImpulseSplitSolution& s_try) {
    s_try.setImpulseStatus(s);
    robot.integrateConfiguration(s.q, d.dq(), step_size, s_try.q);
    s_try.v  = s.v  + step_size * d.dv();
    s_try.dv = s.dv + step_size * d.ddv();
    s_try.f_stack() = s.f_stack() + step_size * d.df();
    s_try.set_f_vector();
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

  static const Eigen::VectorXd& q_prev(const OCPDiscretizer& discretizer, 
                                       const Eigen::VectorXd& q, 
                                       const Solution& s, 
                                       const int time_stage) {
    assert(time_stage >= 0);
    assert(time_stage <= discretizer.N());
    if (time_stage == 0) {
      return q;
    }
    else if (discretizer.isTimeStageBeforeImpulse(time_stage-1)) {
      return s.aux[discretizer.impulseIndexAfterTimeStage(time_stage-1)].q;
    }
    else if (discretizer.isTimeStageBeforeLift(time_stage-1)) {
      return s.lift[discretizer.liftIndexAfterTimeStage(time_stage-1)].q;
    }
    else {
      return s[time_stage-1].q;
    }
  }

  static const Eigen::VectorXd& v_prev(const OCPDiscretizer& discretizer, 
                                       const Eigen::VectorXd& v, 
                                       const Solution& s, 
                                       const int time_stage) {
    assert(time_stage >= 0);
    assert(time_stage <= discretizer.N());
    if (time_stage == 0) {
      return v;
    }
    else if (discretizer.isTimeStageBeforeImpulse(time_stage-1)) {
      return s.aux[discretizer.impulseIndexAfterTimeStage(time_stage-1)].v;
    }
    else if (discretizer.isTimeStageBeforeLift(time_stage-1)) {
      return s.lift[discretizer.liftIndexAfterTimeStage(time_stage-1)].v;
    }
    else {
      return s[time_stage-1].v;
    }
  }

};

} // namespace idocp 

#endif // IDOCP_LINE_SEARCH_HPP_ 