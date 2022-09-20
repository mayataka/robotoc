#include "robotoc/robot/contact_model_info.hpp"

#include <stdexcept>

namespace robotoc {

ContactModelInfo::ContactModelInfo(const std::string& _frame, 
                                   const double _baumgarte_time_step)
  : frame(_frame),
    baumgarte_position_gain(1.0/(_baumgarte_time_step*_baumgarte_time_step)),
    baumgarte_velocity_gain(2.0/(_baumgarte_time_step)) {
  if (_baumgarte_time_step <= 0.0) {
    throw std::out_of_range("[ContactModelInfo] invalid argument: 'baumgarte_time_step' must be positive!");
  }
}


ContactModelInfo::ContactModelInfo(const std::string& _frame, 
                                   const double _baumgarte_position_gain,
                                   const double _baumgarte_velocity_gain)
  : frame(_frame),
    baumgarte_position_gain(_baumgarte_position_gain),
    baumgarte_velocity_gain(_baumgarte_velocity_gain) {
  if (_baumgarte_position_gain <= 0.0) {
    throw std::out_of_range("[ContactModelInfo] invalid argument: 'baumgarte_position_gain' must be positive!");
  }
  if (_baumgarte_velocity_gain <= 0.0) {
    throw std::out_of_range("[ContactModelInfo] invalid argument: 'baumgarte_velocity_gain' must be positive!");
  }
}

} // namespace robotoc