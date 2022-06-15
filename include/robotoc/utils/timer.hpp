#ifndef ROBOTOC_UTILS_TIMER_HPP_
#define ROBOTOC_UTILS_TIMER_HPP_

#include <chrono>

namespace robotoc {

///
/// @class Timer
/// @brief A timer class to take benchmarks. 
///
class Timer {
public:
  ///
  /// @brief Default constructor. 
  ///
  Timer() {
    tick();
    tick_ = tock_;
  }

  ///
  /// @brief Default destructor. 
  ///
  ~Timer() = default;

  ///
  /// @brief Takes a time clock. 
  ///
  void tick() { tick_ = std::chrono::high_resolution_clock::now(); }

  ///
  /// @brief Takes a time clock. 
  ///
  void tock() { tock_ = std::chrono::high_resolution_clock::now(); }

  ///
  
  /// @brief Returns the time duration (seconds) from tick() to tock(). 
  ///
  double s() const {
    const std::chrono::duration<double, std::ratio<1, 1>> timing = tock_ - tick_;
    return timing.count();
  }

  ///
  /// @brief Returns the time duration (milli seconds) from tick() to tock(). 
  ///
  double ms() const {
    const std::chrono::duration<double, std::milli> timing = tock_ - tick_;
    return timing.count();
  }

  ///
  /// @brief Returns the time duration (nano seconds) from tick() to tock(). 
  ///
  double ns() const {
    const std::chrono::duration<double, std::nano> timing = tock_ - tick_;
    return timing.count();
  }

private:
  std::chrono::high_resolution_clock::time_point tick_, tock_;
};

} // namespace robotoc 

#endif // ROBOTOC_UTILS_TIMER_HPP_