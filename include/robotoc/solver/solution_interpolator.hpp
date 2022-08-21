#ifndef ROBOTOC_SOLUTION_INTERPOLATOR_HPP_
#define ROBOTOC_SOLUTION_INTERPOLATOR_HPP_

#include "robotoc/ocp/solution.hpp"
#include "robotoc/hybrid/time_discretization.hpp"


namespace robotoc {

///
/// @class SolutionInterpolator
/// @brief Solution interpolator. 
///
class SolutionInterpolator {
public:
  ///
  /// @brief Constructor. 
  /// @param[in] robot Robot model.
  /// @param[in] N Number of the discretization grids of the horizon except for 
  /// the discrete events. Must be positive.
  /// @param[in] reserved_num_discrete_events Reserved number of discrete events  
  /// on the horizon. 
  ///
  SolutionInterpolator(const Robot& robot, const int N, 
                       const int reserved_num_discrete_events=1);

  ///
  /// @brief Default constructor. 
  ///
  SolutionInterpolator();

  ///
  /// @brief Destructor. 
  ///
  ~SolutionInterpolator();

  ///
  /// @brief Default copy constructor. 
  ///
  SolutionInterpolator(const SolutionInterpolator&) = default;

  ///
  /// @brief Default copy assign operator. 
  ///
  SolutionInterpolator& operator=(const SolutionInterpolator&) = default;

  ///
  /// @brief Default move constructor. 
  ///
  SolutionInterpolator(SolutionInterpolator&&) noexcept = default;

  ///
  /// @brief Default move assign operator. 
  ///
  SolutionInterpolator& operator=(SolutionInterpolator&&) noexcept = default;

  ///
  /// @brief Reserve the discrete-event data. 
  /// @param[in] robot Robot model.
  /// @param[in] reserved_num_discrete_events Reserved number of discrete events  
  /// on the horizon. OCP data for impulse and lift events are constructed 
  /// according to this value. Must be non-negative.
  ///
  void reserve(const Robot& robot, const int reserved_num_discrete_events);

  ///
  /// @return Reserved size of the discrete-event data. 
  ///
  int reservedNumDiscreteEvents() const {
    return stored_solution_.reservedNumDiscreteEvents();
  }

  ///
  /// @brief Stores the current time-discretization and solution. 
  /// @param[in] time_discretization Time discretization. 
  /// @param[out] solution Solution. 
  ///
  void store(const TimeDiscretization& time_discretization,
             const Solution& solution);

  ///
  /// @brief Interpolates the solution. 
  /// @param[in] robot Robot model.
  /// @param[in] time_discretization Time discretization. 
  /// @param[out] solution Solution. 
  ///
  void interpolate(const Robot& robot, 
                   const TimeDiscretization& time_discretization, 
                   Solution& solution) const;

  bool hasStoredSolution() const { return has_stored_solution_; }

private:
  TimeDiscretization stored_time_discretization_;
  Solution stored_solution_;
  bool has_stored_solution_;

  int findStoredImpulseIndexBeforeTime(const double t) const {
    const int N = stored_time_discretization_.N();
    for (int i=1; i<=N; ++i) {
      if (stored_time_discretization_.isTimeStageAfterImpulse(i)) {
        const int impulse_index 
            = stored_time_discretization_.impulseIndexAfterTimeStage(i-1);
        if (t > stored_time_discretization_.gridInfoImpulse(impulse_index).t 
            && t < stored_time_discretization_.gridInfo(i).t) {
          return impulse_index;
        }
      }
    }
    return -1;
  }

  int findStoredLiftIndexBeforeTime(const double t) const {
    const int N = stored_time_discretization_.N();
    for (int i=1; i<=N; ++i) {
      if (stored_time_discretization_.isTimeStageAfterLift(i)) {
        const int lift_index 
            = stored_time_discretization_.liftIndexAfterTimeStage(i-1);
        if (t > stored_time_discretization_.gridInfoLift(lift_index).t 
            && t < stored_time_discretization_.gridInfo(i).t) {
          return lift_index;
        }
      }
    }
    return -1;
  }

  int findStoredGridIndexBeforeTime(const double t) const {
    if (t < stored_time_discretization_.gridInfo(0).t) return -1;
    const int N = stored_time_discretization_.N();
    for (int i=0; i<N; ++i) {
      if (stored_time_discretization_.isTimeStageBeforeImpulse(i)) {
        const int impulse_index 
            = stored_time_discretization_.impulseIndexAfterTimeStage(i);
        if (t < stored_time_discretization_.gridInfoImpulse(impulse_index).t) {
          return i;
        }
      }
      else if (stored_time_discretization_.isTimeStageBeforeLift(i)) {
        const int lift_index 
            = stored_time_discretization_.liftIndexAfterTimeStage(i);
        if (t < stored_time_discretization_.gridInfoLift(lift_index).t) {
          return i;
        }
      } 
      else if (t < stored_time_discretization_.gridInfo(i+1).t) {
        return i;
      }
    }
    return N;
  }

  static void interpolate(const Robot& robot, const SplitSolution& s1, 
                          const SplitSolution& s2, const double alpha, 
                          SplitSolution& s) {
    assert(alpha >= 0.0);
    assert(alpha <= 1.0);
    robot.interpolateConfiguration(s1.q, s2.q, alpha, s.q);
    s.v = (1.0 - alpha) * s1.v + alpha * s2.v;
    s.a = (1.0 - alpha) * s1.a + alpha * s2.a;
    s.u = (1.0 - alpha) * s1.u + alpha * s2.u;
    for (size_t i=0; i<s1.f.size(); ++i) {
      if (s2.isContactActive(i)) 
        s.f[i] = (1.0 - alpha) * s1.f[i] + alpha * s2.f[i];
      else
        s.f[i] = s1.f[i];
    }
    s.lmd  = (1.0 - alpha) * s1.lmd + alpha * s2.lmd;
    s.gmm  = (1.0 - alpha) * s1.gmm + alpha * s2.gmm;
    s.beta = (1.0 - alpha) * s1.beta + alpha * s2.beta;
    for (size_t i=0; i<s1.mu.size(); ++i) {
      if (s2.isContactActive(i)) 
        s.mu[i] = (1.0 - alpha) * s1.mu[i] + alpha * s2.mu[i];
      else
        s.mu[i] = s1.mu[i];
    }
    s.nu_passive = (1.0 - alpha) * s1.nu_passive + alpha * s2.nu_passive;
    s.setContactStatus(s1);
    s.set_f_stack();
    s.set_mu_stack();
    s.setImpulseStatus(s1);
    if (s.hasActiveImpulse()) {
      s.xi_stack() = s1.xi_stack();
    }
  }

  template <typename SplitSolutionType>
  static void interpolatePartial(const Robot& robot, const SplitSolution& s1, 
                                 const SplitSolutionType& s2, const double alpha, 
                                 SplitSolution& s) {
    assert(alpha >= 0.0);
    assert(alpha <= 1.0);
    robot.interpolateConfiguration(s1.q, s2.q, alpha, s.q);
    s.v = (1.0 - alpha) * s1.v + alpha * s2.v;
    s.a = s1.a;
    s.u = s1.u;
    for (size_t i=0; i<s1.f.size(); ++i) {
      s.f[i] = s1.f[i];
    }
    s.lmd  = (1.0 - alpha) * s1.lmd + alpha * s2.lmd;
    s.gmm  = (1.0 - alpha) * s1.gmm + alpha * s2.gmm;
    s.beta = s1.beta;
    for (size_t i=0; i<s1.mu.size(); ++i) {
      s.mu[i] = s1.mu[i];
    }
    s.nu_passive = s1.nu_passive;
    s.setContactStatus(s1);
    s.set_f_stack();
    s.set_mu_stack();
    s.setImpulseStatus(s1);
    if (s.hasActiveImpulse()) {
      s.xi_stack() = s1.xi_stack();
    }
  }

};

} // namespace robotoc

#endif // ROBOTOC_SOLUTION_INTERPOLATOR_HPP_