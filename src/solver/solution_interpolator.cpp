#include "robotoc/solver/solution_interpolator.hpp"

#include <stdexcept>


namespace robotoc {

SolutionInterpolator::SolutionInterpolator(const InterpolationOrder order) 
  : order_(order),
    stored_time_discretization_(),
    stored_solution_(),
    has_stored_solution_(false) {
}


void SolutionInterpolator::setInterpolationOrder(const InterpolationOrder order) {
  order_ = order;
}


void SolutionInterpolator::store(const TimeDiscretization& time_discretization,
                                 const Solution& solution) {
  assert(solution.size() >= time_discretization.size());
  stored_time_discretization_ = time_discretization;
  stored_solution_ = solution;
  has_stored_solution_ = true;
}


void SolutionInterpolator::interpolate(
    const Robot& robot, const TimeDiscretization& time_discretization, 
    Solution& solution) const {
  assert(solution.size() >= time_discretization.size());
  if (!has_stored_solution_) return;

  const int N = time_discretization.size() - 1;
  for (int i=0; i<=N; ++i) {
    const auto& grid = time_discretization[i];
    if (grid.t <= stored_time_discretization_.front().t) {
      solution[i] = stored_solution_[0];
      continue;
    }
    if (grid.t >= stored_time_discretization_.back().t) {
      solution[i] = stored_solution_[stored_time_discretization_.size()-1];
      continue;
    }

    if (grid.type == GridType::Impact) {
      const int stored_grid_index = findStoredGridIndexAtImpactByTime(grid.t);
      if (stored_grid_index >= 0) {
        solution[i] = stored_solution_[stored_grid_index];
        modifyImpactSolution(solution[i]);
        if ((i-2 >= 0) && (stored_grid_index-2 >= 0)) {
          solution[i-2].setSwitchingConstraintDimension(
              stored_solution_[stored_grid_index-2].dims());
          solution[i-2].xi_stack() = stored_solution_[stored_grid_index-2].xi_stack();
        }
      }
      else {
        const int grid_index = findStoredGridIndexBeforeTime(grid.t);
        const double alpha = (grid.t - stored_time_discretization_[grid_index].t) 
                              / stored_time_discretization_[grid_index].dt;
        if (stored_time_discretization_[grid_index+1].type == GridType::Terminal) {
          interpolatePartial(robot, stored_solution_[grid_index], 
                             stored_solution_[grid_index+1], alpha, solution[i]);
          modifyImpactSolution(solution[i]);
        }
        else {
          initEventSolution(robot, stored_solution_[grid_index], 
                            stored_solution_[grid_index+1], alpha, solution[i]);
        }
      }
      continue;
    }

    if (grid.type == GridType::Lift) {
      const int stored_grid_index = findStoredGridIndexAtLiftByTime(grid.t);
      if (stored_grid_index >= 0) {
        solution[i] = stored_solution_[stored_grid_index];
      }
      else {
        const int grid_index = findStoredGridIndexBeforeTime(grid.t);
        const double alpha = (grid.t - stored_time_discretization_[grid_index].t) 
                              / stored_time_discretization_[grid_index].dt;
        if (stored_time_discretization_[grid_index+1].type == GridType::Terminal) {
          interpolatePartial(robot, stored_solution_[grid_index], 
                             stored_solution_[grid_index+1], alpha, solution[i]);
        }
        else {
          initEventSolution(robot, stored_solution_[grid_index], 
                            stored_solution_[grid_index+1], alpha, solution[i]);
        }
      }
      continue;
    }

    const int grid_index = findStoredGridIndexBeforeTime(grid.t);
    const double alpha = (grid.t - stored_time_discretization_[grid_index].t) 
                          / stored_time_discretization_[grid_index].dt;
    if (order_ == InterpolationOrder::Zero) {
      solution[i] = stored_solution_[grid_index];
      continue;
    }
    if (stored_time_discretization_[grid_index+1].type != GridType::Intermediate) {
      interpolatePartial(robot, stored_solution_[grid_index],
                         stored_solution_[grid_index+1], alpha, solution[i]);
      continue;
    }
    interpolate(robot, stored_solution_[grid_index],
                stored_solution_[grid_index+1], alpha, solution[i]);
  }
  modifyTerminalSolution(solution[N]);
}


void SolutionInterpolator::interpolate(const Robot& robot, 
                                       const SplitSolution& s1, 
                                       const SplitSolution& s2, 
                                       const double alpha, SplitSolution& s) {
  assert(alpha >= 0.0);
  assert(alpha <= 1.0);
  robot.interpolateConfiguration(s1.q, s2.q, alpha, s.q);
  s.v = (1.0 - alpha) * s1.v + alpha * s2.v;
  s.u = (1.0 - alpha) * s1.u + alpha * s2.u;
  s.a = (1.0 - alpha) * s1.a + alpha * s2.a;
  s.dv.setZero();
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
  s.set_f_stack();
  s.set_mu_stack();
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
  s.u = s1.u;
  s.a = s1.a;
  s.dv.setZero();
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
  s.set_f_stack();
  s.set_mu_stack();
}


void SolutionInterpolator::initEventSolution(const Robot& robot, 
                                             const SplitSolution& s1, 
                                             const SplitSolution& s2, 
                                             const double alpha, 
                                             SplitSolution& s) {
  assert(alpha >= 0.0);
  assert(alpha <= 1.0);
  robot.interpolateConfiguration(s1.q, s2.q, alpha, s.q);
  s.v = (1.0 - alpha) * s1.v + alpha * s2.v;
  s.u = s2.u;
  s.a = (1.0 - alpha) * s1.a + alpha * s2.a;
  s.dv.setZero();
  for (size_t i=0; i<s2.f.size(); ++i) {
    s.f[i] = s2.f[i];
  }
  s.lmd  = (1.0 - alpha) * s1.lmd + alpha * s2.lmd;
  s.gmm  = (1.0 - alpha) * s1.gmm + alpha * s2.gmm;
  s.beta = s2.beta;
  for (size_t i=0; i<s2.mu.size(); ++i) {
    s.mu[i] = s2.mu[i];
  }
  s.nu_passive = s2.nu_passive;
  s.set_f_stack();
  s.set_mu_stack();
}


void SolutionInterpolator::modifyImpactSolution(SplitSolution& s) {
  s.u.setZero();
  s.a.setZero();
  s.nu_passive.setZero();
}


void SolutionInterpolator::modifyTerminalSolution(SplitSolution& s) {
  s.u.setZero();
  s.a.setZero();
  for (auto& e : s.f) {
    e.setZero();
  }
  s.beta.setZero();
  for (auto& e : s.mu) {
    e.setZero();
  }
  s.nu_passive.setZero();
  s.set_f_stack();
  s.set_mu_stack();
}

} // namespace robotoc