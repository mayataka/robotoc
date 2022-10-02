#include "solution_factory.hpp"

#include "robotoc/ocp/time_discretization.hpp"


namespace robotoc {
namespace testhelper {

Solution CreateSolution(const Robot& robot, const int N, 
                        const int max_num_impulse) {
  Solution s(robot, N, max_num_impulse);
  for (int i=0; i<=N; ++i) {
    s[i].setRandom(robot);
  }
  return s;
}


Solution CreateSolution(const Robot& robot,  
                        const std::shared_ptr<ContactSequence>& contact_sequence, 
                        const double T, const int N, 
                        const int max_num_impulse, const double t) {
  if (robot.maxNumContacts() == 0) {
    return CreateSolution(robot, N, max_num_impulse);
  }
  else {
    TimeDiscretization time_discretization(T, N, max_num_impulse);
    time_discretization.discretize(contact_sequence, t);
    Solution s(robot, time_discretization.N_grids(), max_num_impulse);
    for (int i=0; i<time_discretization.N_grids(); ++i) {
      const auto& grid = time_discretization.grid(i);
      if (grid.type == GridType::Impulse) {
        s[i].setContactStatus(contact_sequence->impulseStatus(grid.impulse_index));
        s[i].setRandom(robot);
      }
      else {
        s[i].setContactStatus(contact_sequence->contactStatus(grid.contact_phase));
        if (grid.switching_constraint) {
          s[i].setSwitchingConstraintDimension(contact_sequence->impulseStatus(grid.impulse_index+1).dimf());
        }
        else {
          s[i].setSwitchingConstraintDimension(0);
        }
        s[i].setRandom(robot);
      }
    }
    s[time_discretization.N_grids()].setRandom(robot);
    return s;
  }
}

} // namespace testhelper
} // namespace robotoc