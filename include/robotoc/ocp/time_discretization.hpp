#ifndef ROBOTOC_TIME_DISCRETIZATION_HPP_ 
#define ROBOTOC_TIME_DISCRETIZATION_HPP_

#include <vector>
#include <memory>
#include <limits>
#include <cmath>
#include <iostream>

#include "robotoc/planner/contact_sequence.hpp"
#include "robotoc/ocp/grid_info.hpp"


namespace robotoc {

///
/// @class TimeDiscretization
/// @brief Time discretization of the hybrid optimal control problem.
///
class TimeDiscretization {
public:
  ///
  /// @brief Constructor. 
  /// @param[in] T Length of the horizon. Must be positive.
  /// @param[in] N Number of the discretization grids of the horizon except for 
  /// the discrete events. Must be positive.
  /// @param[in] reserved_num_discrete_events Reserved size of each discrete 
  /// events (impulse and lift) to avoid dynamic memory allocation. Must be 
  /// non-negative. Default is 0.
  ///
  TimeDiscretization(const double T, const int N, 
                     const int reserved_num_discrete_events=0);

  ///
  /// @brief Default constructor. 
  ///
  TimeDiscretization();

  ///
  /// @brief Default destructor. 
  ///
  ~TimeDiscretization() = default;

  ///
  /// @brief Default copy constructor. 
  ///
  TimeDiscretization(const TimeDiscretization&) = default;

  ///
  /// @brief Default copy assign operator. 
  ///
  TimeDiscretization& operator=(const TimeDiscretization&) = default;

  ///
  /// @brief Default move constructor. 
  ///
  TimeDiscretization(TimeDiscretization&&) noexcept = default;

  ///
  /// @brief Default move assign operator. 
  ///
  TimeDiscretization& operator=(TimeDiscretization&&) noexcept = default;

  ///
  /// @return Number of the time stages on the horizon. 
  ///
  int N() const {
    return N_;
  }

  ///
  /// @brief Returns the number of grids. 
  /// @return The number of grids..
  ///
  inline int size() const {
    return num_grids_ + 1;
  }

  ///
  /// @brief Returns the grid info of the specified stage. 
  /// @param[in] i Stage of interest. 
  /// @return const reference to the grid info of the stage of interest.
  ///
  inline const GridInfo& grid(const int i) const {
    assert(i >= 0);
    assert(i < size());
    return grid_[i];
  }

  ///
  /// @brief Returns the grid info of the specified stage. 
  /// @param[in] i Stage of interest. 
  /// @return const reference to the grid info of the stage of interest.
  ///
  inline const GridInfo& operator[] (const int i) const {
    assert(i >= 0);
    assert(i < size());
    return grid_[i];
  }

  ///
  /// @brief Returns the grid info of the initial stage. 
  /// @return const reference to the grid info of the intial stage.
  ///
  inline const GridInfo& front() const {
    return grid_[0];
  }

  ///
  /// @brief Returns the grid info of the terminal stage. 
  /// @return const reference to the grid info of the terminal stage.
  ///
  inline const GridInfo& back() const {
    return grid_[size()-1];
  }

  ///
  /// @brief Gets the maximum time step of the discretization.
  /// @return The maximum time step of the discretization.
  ///
  double maxTimeStep() const {
    double max_dt = grid_[0].dt;
    for (int i=0; i<num_grids_; ++i) {
      max_dt = std::max(grid_[i].dt, max_dt);
    }
    return max_dt;
  }

  ///
  /// @brief Reserve the discrete-event data. 
  /// @param[in] reserved_num_discrete_events Reserved size of discrete events  
  /// on the horizon. Must be non-negative.
  ///
  void reserve(const int reserved_num_discrete_events);

  ///
  /// @return Reserved size of the discrete-event data. 
  ///
  int reservedNumDiscreteEvents() const;

  ///
  /// @brief Discretizes the finite horizon taking into account the discrete 
  /// events.
  /// @param[in] contact_sequence Shared ptr to the contact sequence.
  /// @param[in] t Initial time of the horizon.
  ///
  void discretize(const std::shared_ptr<ContactSequence>& contact_sequence, const double t);

  ///
  /// @brief Discretizes the finite horizon taking into account the discrete 
  /// events.
  /// @param[in] contact_sequence Shared ptr to the contact sequence.
  /// @param[in] t Initial time of the horizon.
  ///
  void correctTimeSteps(const std::shared_ptr<ContactSequence>& contact_sequence, const double t);

  ///
  /// @brief Displays the time discretization onto a ostream.
  ///
  void disp(std::ostream& os) const;

  friend std::ostream& operator<<(std::ostream& os, 
                                  const TimeDiscretization& discretization);

private:
  double T_, max_dt_, eps_;
  int N_, num_grids_, reserved_num_discrete_events_;
  std::vector<GridInfo> grid_;
  std::vector<bool> sto_event_, sto_phase_;
};

} // namespace robotoc

#endif // ROBOTOC_TIME_DISCRETIZATION_HPP_ 