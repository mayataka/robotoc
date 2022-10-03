#include "direction_factory.hpp"


namespace robotoc {
namespace testhelper {

Direction CreateDirection(const Robot& robot, const int N) {
  Direction d(N+1, SplitDirection(robot));
  for (int i=0; i<=N; ++i) {
    d[i].setRandom();
  }
  return d;
}


Direction CreateDirection(const Robot& robot, 
                          const std::shared_ptr<ContactSequence>& contact_sequence,
                          const TimeDiscretization& time_discretization) {
  if (robot.maxNumContacts() == 0) {
    return CreateDirection(robot, time_discretization.N());
  }
  else {
    Direction d(time_discretization.size(), SplitDirection(robot));
    for (int i=0; i<time_discretization.size()-1; ++i) {
      const auto& grid = time_discretization[i];
      if (grid.type == GridType::Impulse) {
        d[i].setContactDimension(contact_sequence->impulseStatus(grid.impulse_index).dimf());
        d[i].setRandom();
      }
      else {
        d[i].setContactDimension(contact_sequence->contactStatus(grid.phase).dimf());
        if (grid.switching_constraint) {
          d[i].setSwitchingConstraintDimension(contact_sequence->impulseStatus(grid.impulse_index+1).dimf());
        }
        else {
          d[i].setSwitchingConstraintDimension(0);
        }
        d[i].setRandom();
      }
    }
    d[time_discretization.size()-1].setRandom();
    return d;
  }
}

} // namespace testhelper
} // namespace robotoc