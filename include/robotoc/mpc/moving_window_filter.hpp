#ifndef ROBOTOC_MOVING_WINDOW_FILTER_HPP_
#define ROBOTOC_MOVING_WINDOW_FILTER_HPP_

#include <stdexcept>
#include <iostream>
#include <cassert>

#include "robotoc/utils/aligned_deque.hpp"

namespace robotoc {

///
/// @class MovingWindowFilter
/// @brief Moving window filter for foot step planning. 
///
template <int dim>
class MovingWindowFilter {
public:
  using Vector = Eigen::Matrix<double, dim, 1>;

  ///
  /// @brief Constructs the filter.
  /// @param[in] time_length Time length (seconds) of the filter. 
  /// @param[in] min_sampling_period Minimum sampling period. Must be 
  /// non-negative. Default is zero.
  ///
  MovingWindowFilter(const double time_length, 
                     const double min_sampling_period=0.0)
    : time_length_(time_length),
      min_sampling_period_(min_sampling_period),
      last_sampling_time_(0.0),
      time_(),
      data_(),
      average_(Vector::Zero()) {
    try {
      if (time_length <= 0.0) {
        throw std::out_of_range("invalid argument: time_length must be positive!");
      }
      if (min_sampling_period < 0) {
        throw std::out_of_range("invalid argument: min_sampling_period must be non-negative!");
      }
    }
    catch(const std::exception& e) {
      std::cerr << e.what() << '\n';
      std::exit(EXIT_FAILURE);
    }
  }

  ///
  /// @brief Default constructor. 
  ///
  MovingWindowFilter() 
    : time_length_(0),
      min_sampling_period_(0.0),
      last_sampling_time_(0.0),
      time_(),
      data_(),
      average_(Vector::Zero()) {
  }

  ///
  /// @brief Destructor. 
  ///
  ~MovingWindowFilter() {}

  ///
  /// @brief Default copy constructor. 
  ///
  MovingWindowFilter(const MovingWindowFilter&) = default;

  ///
  /// @brief Default copy assign operator. 
  ///
  MovingWindowFilter& operator=(const MovingWindowFilter&) = default;

  ///
  /// @brief Default move constructor. 
  ///
  MovingWindowFilter(MovingWindowFilter&&) noexcept = default;

  ///
  /// @brief Default move assign operator. 
  ///
  MovingWindowFilter& operator=(MovingWindowFilter&&) noexcept = default;

  ///
  /// @brief Set parameters of the filter.
  /// @param[in] time_length Time length (seconds) of the filter. 
  /// @param[in] min_sampling_period Minimum sampling period. Must be 
  /// non-negative. Default is zero.
  ///
  void setParameters(const double time_length, 
                     const double min_sampling_period=0.0) {
    try {
      if (time_length <= 0.0) {
        throw std::out_of_range("invalid argument: time_length must be positive!");
      }
      if (min_sampling_period < 0) {
        throw std::out_of_range("invalid argument: min_sampling_period must be non-negative!");
      }
    }
    catch(const std::exception& e) {
      std::cerr << e.what() << '\n';
      std::exit(EXIT_FAILURE);
    }
    time_length_ = time_length;
    min_sampling_period_ = min_sampling_period;
  }

  ///
  /// @brief Clear the filter.
  ///
  void clear() {
    last_sampling_time_ = 0.0;
    time_.clear();
    data_.clear();
    average_.setZero();
  }

  ///
  /// @brief Push back a data.
  /// @param[in] t Time of the data sampling.
  /// @param[in] data The data.
  ///
  void push_back(const double t, const Vector& data) {
    if (time_.empty()) {
      time_.push_back(t);
      data_.push_back(data);
      last_sampling_time_ = t;
    }
    else {
      if (t - last_sampling_time_ >= min_sampling_period_) {
        time_.push_back(t);
        data_.push_back(data);
        while (t - time_.front() > time_length_) {
          time_.pop_front();
          data_.pop_front();
        }
        average_.setZero();
        for (const auto& e : data_) {
          average_.noalias() += e;
        }
        average_.array() /= data_.size();
        last_sampling_time_ = t;
      }
    }
  }

  ///
  /// @brief Gets the current data size.
  ///
  int size() const {
    return data_.size();
  }

  ///
  /// @brief Gets the average.
  ///
  const Vector& average() const {
    return average_;
  }

  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

private:
  double time_length_, min_sampling_period_, last_sampling_time_;
  std::deque<double> time_;
  aligned_deque<Vector> data_;
  Vector average_;

};

} // namespace robotoc 

#endif // ROBOTOC_MOVING_WINDOW_FILTER_HPP_