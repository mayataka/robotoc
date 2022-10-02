#include "robotoc/solver/solution_interpolator.hpp"


namespace robotoc {

SolutionInterpolator::SolutionInterpolator(const Robot& robot, const int N, 
                                           const int reserved_num_discrete_events) 
  : stored_time_discretization_(),
    stored_solution_(N+1+3*reserved_num_discrete_events, SplitSolution(robot)),
    has_stored_solution_(false) {
}


SolutionInterpolator::SolutionInterpolator()
  : stored_time_discretization_(),
    stored_solution_(),
    has_stored_solution_(false) {
}


void SolutionInterpolator::reserve(const Robot& robot, 
                                   const TimeDiscretization& time_discretization) {
  while (stored_solution_.size() < time_discretization.N_grids()+1) {
    stored_solution_.push_back(stored_solution_.back());
  }
}


void SolutionInterpolator::store(const TimeDiscretization& time_discretization,
                                 const Solution& solution) {
  stored_time_discretization_ = time_discretization;
  stored_solution_ = solution;
  has_stored_solution_ = true;
}


void SolutionInterpolator::interpolate(
    const Robot& robot, const TimeDiscretization& time_discretization, 
    Solution& solution) const {
  if (!has_stored_solution_) return;

  const int N = time_discretization.N_grids();
  int next_impulse_index = time_discretization.grid(0).impulse_index + 1;
  int next_lift_index = time_discretization.grid(0).lift_index + 1;
  for (int i=0; i<=N; ++i) {
    const auto& grid = time_discretization.grid(i);
    const auto& grid_next = time_discretization.grid(i+1);
    if (grid.t <= stored_time_discretization_.grid(0).t) {
      solution[i] = stored_solution_[0];
      continue;
    }
    if (grid.t >= stored_time_discretization_.grid(N).t) {
      solution[i] = stored_solution_[N];
      continue;
    }
    if (grid.type == GridType::Impulse) {
      const int grid_index = findStoredGridIndexAtImpulse(grid.t);
      assert(grid_index >= 0);
      solution[i] = stored_solution_[grid_index];
      continue;
    }
    else if (grid.type == GridType::Lift) {
      const int grid_index = findStoredGridIndexAtLift(grid.t);
      assert(grid_index >= 0);
      solution[i] = stored_solution_[grid_index];
      continue;
    }
    else {
      const int grid_index = findStoredGridIndexBeforeTime(grid.t);
      const double alpha = (grid.t - stored_time_discretization_.grid(grid_index).t) 
                            / stored_time_discretization_.grid(grid_index).dt;
      if ((stored_time_discretization_.grid(grid_index+1).type == GridType::Impulse)
          || (stored_time_discretization_.grid(grid_index+1).type == GridType::Lift)) {
        interpolatePartial(robot, stored_solution_[grid_index],
                           stored_solution_[grid_index+1], alpha, solution[i]);
      }
      else {
        interpolate(robot, stored_solution_[grid_index],
                    stored_solution_[grid_index+1], alpha, solution[i]);
      }
    }
  }
}


void SolutionInterpolator::interpolate(const Robot& robot, 
                                       const SplitSolution& s1, 
                                       const SplitSolution& s2, 
                                       const double alpha, SplitSolution& s) {
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
  s.setSwitchingConstraintDimension(s1.dims());
  if (s.dims() > 0) {
    s.xi_stack() = s1.xi_stack();
  }
}


void SolutionInterpolator::interpolatePartial(const Robot& robot, 
                                              const SplitSolution& s1, 
                                              const SplitSolution& s2, 
                                              const double alpha, 
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
  s.setSwitchingConstraintDimension(s1.dims());
  if (s.dims() > 0) {
    s.xi_stack() = s1.xi_stack();
  }
}

} // namespace robotoc