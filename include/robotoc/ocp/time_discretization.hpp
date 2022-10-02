#ifndef ROBOTOC_TIME_DISCRETIZATION_HPP_ 
#define ROBOTOC_TIME_DISCRETIZATION_HPP_

#include <vector>
#include <memory>
#include <limits>
#include <cmath>
#include <iostream>

#include "robotoc/planner/contact_sequence.hpp"
#include "robotoc/ocp/discretization_method.hpp"
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
  /// @brief Sets the discretization method of the optimal contro problem. 
  /// @param[in] discretization_method The discretization method.
  ///
  void setDiscretizationMethod(const DiscretizationMethod discretization_method);

  ///
  /// @return Number of the time stages on the horizon. 
  ///
  int N() const {
    return N_;
  }

  int N_grids() const {
    return N_grids_;
  }

  ///
  /// @brief Returns the grid info of the specified stage. 
  /// @param[in] i Stage of interest. 
  /// @return const reference to the grid info of the stage of interest.
  ///
  inline const GridInfo& gridInfo(const int i) const {
    assert(i >= 0);
    assert(i <= N_grids());
    return grid_[i];
  }

  ///
  /// @brief Returns the grid info of the specified stage. 
  /// @param[in] i Stage of interest. 
  /// @return const reference to the grid info of the stage of interest.
  ///
  inline const GridInfo& grid(const int i) const {
    assert(i >= 0);
    assert(i <= N_grids());
    return grid_[i];
  }

  ///
  /// @brief Returns the grid info of the specified stage. 
  /// @param[in] i Stage of interest. 
  /// @return const reference to the grid info of the stage of interest.
  ///
  inline const GridInfo& operator[] (const int i) const {
    assert(i >= 0);
    assert(i <= N_grids());
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
    return grid_[N_grids()];
  }

  ///
  /// @brief Returns the number of grids. 
  /// @return The number of grids..
  ///
  inline int size() const {
    return N_grids() + 1;
  }

  ///
  /// @brief Returns the current discretization method. 
  /// @return The current discretization method.
  ///
  DiscretizationMethod discretizationMethod() const;

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

  void discretizeGrid(const std::shared_ptr<ContactSequence>& contact_sequence, const double t);

  void discretizePhase(const std::shared_ptr<ContactSequence>& contact_sequence, const double t);

  ///
  /// @brief Discretizes the finite horizon taking into account the discrete 
  /// events.
  /// @param[in] contact_sequence Shared ptr to the contact sequence.
  /// @param[in] t Initial time of the horizon.
  /// @note If the discretization method is DiscretizationMethod::GridBased, 
  /// this funtion can change the structure of the discretization, i.e., 
  /// the number of grids on each contact phase. If the discretization method is 
  /// DiscretizationMethod::PhaseBased, this function keeps the structure of the 
  /// discretization. In the latter case, meshRefinement() is needed to chagne 
  /// the discretization structure.
  ///
  void discretize(const std::shared_ptr<ContactSequence>& contact_sequence, 
                  const double t) {
    discretizeGrid(contact_sequence, t);
  }

  ///
  /// @brief Applies the mesh refinement and changes the structure of the 
  /// discretization, i.e., the number of grids on each contact phase for 
  /// the discretization method DiscretizationMethod::PhaseBased.
  /// Specifically, this function reduces the numbers of the grids
  /// from the phases where the solution is relatively accurate and increases
  /// them to the phases where the solution is relatively inaccurate while 
  /// keeping the total number of the discretization grids including the 
  /// discrete events. 
  /// @param[in] contact_sequence Contact sequence.
  /// @param[in] t Initial time of the horizon.
  /// @note This function do nothing if the discretization method is 
  /// DiscretizationMethod::GridBased.
  ///
  void meshRefinement(const std::shared_ptr<ContactSequence>& contact_sequence, 
                      const double t) {
    discretizeGrid(contact_sequence, t);
    discretizePhase(contact_sequence, t);
  }

  inline const std::vector<GridInfo>& getGrid() const {
    return grid_;
  }

  double dt_max() const {
    double max_dt = grid_[0].dt;
    for (int i=0; i<N_grids(); ++i) {
      max_dt = std::max(grid_[i].dt, max_dt);
    }
    return max_dt;
  }

  ///
  /// @brief Displays the discretization of the hybrid optimal control problem 
  /// onto a ostream.
  ///
  void disp(std::ostream& os) const;

  friend std::ostream& operator<<(std::ostream& os, 
                                  const TimeDiscretization& discretization);

private:
  double T_, max_dt_, eps_;
  int N_, N_grids_, reserved_num_discrete_events_;
  std::vector<GridInfo> grid_;
  std::vector<bool> sto_event_, sto_phase_;
  DiscretizationMethod discretization_method_;
};

} // namespace robotoc

#endif // ROBOTOC_TIME_DISCRETIZATION_HPP_ 