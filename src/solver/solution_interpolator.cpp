#include "robotoc/solver/solution_interpolator.hpp"


namespace robotoc {

SolutionInterpolator::SolutionInterpolator(const Robot& robot, const int N, 
                                           const int reserved_num_discrete_events) 
  : stored_time_discretization_(),
    stored_solution_(robot, N, reserved_num_discrete_events) {
}


SolutionInterpolator::SolutionInterpolator()
  : stored_time_discretization_(),
    stored_solution_() {
}


SolutionInterpolator::~SolutionInterpolator() {
}


void SolutionInterpolator::reserve(const Robot& robot, 
                                   const int reserved_num_discrete_events) {
  stored_solution_.reserve(robot, reserved_num_discrete_events);
}


void SolutionInterpolator::storeTimeDiscretization(
    const TimeDiscretization& time_discretization) {
  stored_time_discretization_ = time_discretization;
}


void SolutionInterpolator::interpolate(
    const Robot& robot, const TimeDiscretization& time_discretization, 
    Solution& solution) {
  stored_solution_ = solution;
  const int N = time_discretization.N();
  for (int i=0; i<N; ++i) {
    const double t = time_discretization.gridInfo(i).t;
    if (t <= stored_time_discretization_.gridInfo(0).t) {
      solution[i] = stored_solution_[0];
      continue;
    }
    if (t >= stored_time_discretization_.gridInfo(N).t) {
      solution[i] = stored_solution_[N];
      continue;
    }
    const int impulse_index = findStoredImpulseIndexBeforeTime(t);
    if (impulse_index >= 0) {
      const int grid_index = stored_time_discretization_.timeStageAfterImpulse(impulse_index);
      const double alpha = (t - stored_time_discretization_.gridInfoImpulse(impulse_index).t) 
                            / stored_time_discretization_.gridInfoImpulse(impulse_index).dt;
      interpolate(robot, stored_solution_.aux[impulse_index],
                  stored_solution_[grid_index], alpha, solution[i]);
    }
    const int lift_index = findStoredLiftIndexBeforeTime(t);
    if (lift_index >= 0) {
      const int grid_index = stored_time_discretization_.timeStageAfterLift(lift_index);
      const double alpha = (t - stored_time_discretization_.gridInfoLift(lift_index).t) 
                            / stored_time_discretization_.gridInfoLift(lift_index).dt;
      interpolate(robot, stored_solution_.lift[lift_index],
                  stored_solution_[grid_index], alpha, solution[i]);
    }
    const int grid_index = findStoredGridIndexBeforeTime(t);
    if (stored_time_discretization_.isTimeStageBeforeImpulse(grid_index)
        || stored_time_discretization_.isTimeStageBeforeLift(grid_index)) {
      solution[i] = stored_solution_[grid_index];
    }
    else {
      const double alpha = (t - stored_time_discretization_.gridInfo(grid_index).t) 
                            / stored_time_discretization_.gridInfo(grid_index).dt;
      interpolate(robot, stored_solution_[grid_index],
                  stored_solution_[grid_index+1], alpha, solution[i]);
    }
  }
}

} // namespace robotoc