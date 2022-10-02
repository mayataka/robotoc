#include "direction_factory.hpp"

#include "robotoc/ocp/time_discretization.hpp"


namespace robotoc {
namespace testhelper {

Direction CreateDirection(const Robot& robot, const int N, 
                          const int max_num_impulse) {
  Direction d(robot, N, max_num_impulse);
  for (int i=0; i<=N; ++i) {
    d[i].setRandom();
  }
  return d;
}


Direction CreateDirection(const Robot& robot, 
                          const std::shared_ptr<ContactSequence>& contact_sequence, 
                          const double T, const int N, 
                          const int max_num_impulse, const double t) {
  if (robot.maxNumContacts() == 0) {
    return CreateDirection(robot, N, max_num_impulse);
  }
  else {
    TimeDiscretization time_discretization(T, N, max_num_impulse);
    time_discretization.discretize(contact_sequence, t);
    Direction d(robot, time_discretization.N_grids(), max_num_impulse);
    for (int i=0; i<time_discretization.N_grids(); ++i) {
      const auto& grid = time_discretization.grid(i);
      if (grid.type == GridType::Impulse) {
        d[i].setContactDimension(contact_sequence->impulseStatus(grid.impulse_index).dimf());
        d[i].setRandom();
      }
      else {
        d[i].setContactDimension(contact_sequence->contactStatus(grid.contact_phase).dimf());
        if (grid.switching_constraint) {
          d[i].setSwitchingConstraintDimension(contact_sequence->impulseStatus(grid.impulse_index+1).dimf());
        }
        else {
          d[i].setSwitchingConstraintDimension(0);
        }
        d[i].setRandom();
      }
    }
    d[time_discretization.N_grids()].setRandom();
    return d;
  }
}

} // namespace testhelper
} // namespace robotoc