#ifndef ROBOTOC_SOLUTION_INTERPOLATOR_HPP_
#define ROBOTOC_SOLUTION_INTERPOLATOR_HPP_

#include "robotoc/core/solution.hpp"
#include "robotoc/ocp/time_discretization.hpp"
#include "robotoc/utils/numerics.hpp"


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
  /// @brief Default destructor. 
  ///
  ~SolutionInterpolator() = default;

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
  /// @brief Stores the current time discretization and solution. 
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
  void interpolateEventBased(const Robot& robot, 
                             const TimeDiscretization& time_discretization, 
                             Solution& solution) const;

  ///
  /// @brief Interpolates the solution. 
  /// @param[in] robot Robot model.
  /// @param[in] time_discretization Time discretization. 
  /// @param[out] solution Solution. 
  ///
  void interpolateTimeBased(const Robot& robot, 
                            const TimeDiscretization& time_discretization, 
                            Solution& solution) const;

  ///
  /// @brief Check if this has a stored solution. 
  /// @return true if this has a stored solution. 
  ///
  bool hasStoredSolution() const { return has_stored_solution_; }

private:
  TimeDiscretization stored_time_discretization_;
  Solution stored_solution_;
  bool has_stored_solution_;

  int findStoredGridIndexAtImpulseByTime(const double t) const {
    const int N = stored_time_discretization_.size() - 1;
    for (int i=1; i<N; ++i) {
      if ((stored_time_discretization_[i].type == GridType::Impulse)
            && (numerics::isApprox(t, stored_time_discretization_[i].t))) {
        return i;
      }
    }
    return -1;
  }

  int findStoredGridIndexAtLiftByTime(const double t) const {
    const int N = stored_time_discretization_.size() - 1;
    for (int i=1; i<N; ++i) {
      if ((stored_time_discretization_[i].type == GridType::Lift)
            && (numerics::isApprox(t, stored_time_discretization_[i].t))) {
        return i;
      }
    }
    return -1;
  }

  int findStoredGridIndexAtImpulseByIndex(const int impulse_index) const {
    const int N = stored_time_discretization_.size() - 1;
    for (int i=1; i<N; ++i) {
      if ((stored_time_discretization_[i].type == GridType::Impulse)
            && (impulse_index == stored_time_discretization_[i].impulse_index)) {
        return i;
      }
    }
    return -1;
  }

  int findStoredGridIndexAtLiftByIndex(const int lift_index) const {
    const int N = stored_time_discretization_.size() - 1;
    for (int i=1; i<N; ++i) {
      if ((stored_time_discretization_[i].type == GridType::Lift)
            && (lift_index == stored_time_discretization_[i].lift_index)) {
        return i;
      }
    }
    return -1;
  }

  int findStoredGridIndexBeforeTime(const double t) const {
    if (t < stored_time_discretization_[0].t) return -1;
    const int N = stored_time_discretization_.size() - 1;
    for (int i=0; i<N; ++i) {
      if (t < stored_time_discretization_[i+1].t) {
        return i;
      }
    }
    return N;
  }

  static void interpolate(const Robot& robot, const SplitSolution& s1, 
                          const SplitSolution& s2, const double alpha, 
                          SplitSolution& s);

  static void interpolatePartial(const Robot& robot, const SplitSolution& s1, 
                                 const SplitSolution& s2, const double alpha, 
                                 SplitSolution& s);
};

} // namespace robotoc

#endif // ROBOTOC_SOLUTION_INTERPOLATOR_HPP_