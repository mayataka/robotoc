#include "contact_sequence_factory.hpp"

#include "idocp/hybrid/discrete_event.hpp"


namespace idocp {
namespace testhelper {

ContactSequence CreateContactSequence(const Robot& robot, const int N, 
                                      const int max_num_impulse,
                                      const double t0,
                                      const double event_period) {
  if (robot.maxPointContacts() > 0) {
    std::vector<DiscreteEvent> discrete_events;
    std::vector<double> event_times;
    ContactStatus pre_contact_status = robot.createContactStatus();
    pre_contact_status.setRandom();
    ContactSequence contact_sequence(robot, max_num_impulse);
    contact_sequence.setContactStatusUniformly(pre_contact_status);
    ContactStatus post_contact_status = pre_contact_status;
    std::random_device rnd;
    for (int i=0; i<max_num_impulse; ++i) {
      DiscreteEvent tmp(robot.maxPointContacts());
      tmp.setDiscreteEvent(pre_contact_status, post_contact_status);
      while (!tmp.existDiscreteEvent()) {
        post_contact_status.setRandom();
        tmp.setDiscreteEvent(pre_contact_status, post_contact_status);
      }
      discrete_events.push_back(tmp);
      const double event_time = t0 + i * event_period + 0.1 * event_period * std::abs(Eigen::VectorXd::Random(1)[0]);
      event_times.push_back(event_time);
      pre_contact_status = post_contact_status;
    }
    for (int i=0; i<max_num_impulse; ++i) {
      contact_sequence.push_back(discrete_events[i], event_times[i]);
    }
    return contact_sequence;
  }
  else {
    return ContactSequence(robot, N);
  }
}
  
} // namespace testhelper
} // namespace idocp