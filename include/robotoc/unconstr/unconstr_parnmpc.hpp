#ifndef ROBOTOC_UNCONSTR_PARNMPC_HPP_
#define ROBOTOC_UNCONSTR_PARNMPC_HPP_

#include <vector>
#include <memory>
#include <stdexcept>
#include <iostream>

#include "robotoc/robot/robot.hpp"
#include "robotoc/unconstr/split_unconstr_parnmpc.hpp"
#include "robotoc/unconstr/terminal_unconstr_parnmpc.hpp"
#include "robotoc/cost/cost_function.hpp"
#include "robotoc/constraints/constraints.hpp"


namespace robotoc {

///
/// @class UnconstrParNMPC
/// @brief An optimal control problem of unconstrained rigid-body systems for 
/// ParNMPC algorithm.
///
class UnconstrParNMPC {
public:
  ///
  /// @brief Constructor. 
  /// @param[in] robot Robot model. 
  /// @param[in] cost Shared ptr to the cost function.
  /// @param[in] constraints Shared ptr to the constraints.
  /// @param[in] T Length of the horizon.
  /// @param[in] N Number of the discretization grids of the horizon.
  ///
  UnconstrParNMPC(const Robot& robot, const std::shared_ptr<CostFunction>& cost, 
                  const std::shared_ptr<Constraints>& constraints, 
                  const double T, const int N) 
    : data(N-1, SplitUnconstrParNMPC(robot, cost, constraints)), 
      terminal(TerminalUnconstrParNMPC(robot, cost, constraints)),
      robot_(robot),
      cost_(cost),
      constraints_(constraints),
      T_(T),
      N_(N) {
    try {
      if (T <= 0) {
        throw std::out_of_range("invalid value: T must be positive!");
      }
      if (N <= 0) {
        throw std::out_of_range("invalid value: N must be positive!");
      }
    }
    catch(const std::exception& e) {
      std::cerr << e.what() << '\n';
      std::exit(EXIT_FAILURE);
    }
  }

  ///
  /// @brief Default Constructor.
  ///
  UnconstrParNMPC() 
    : data(), 
      terminal(),
      robot_(),
      cost_(),
      constraints_(),
      T_(0),
      N_(0) {
  }

  ///
  /// @brief Destructor.
  ///
  ~UnconstrParNMPC() {
  }

  ///
  /// @brief Default copy constructor. 
  ///
  UnconstrParNMPC(const UnconstrParNMPC&) = default;

  ///
  /// @brief Default copy assign operator. 
  ///
  UnconstrParNMPC& operator=(const UnconstrParNMPC&) = default;

  ///
  /// @brief Default move constructor. 
  ///
  UnconstrParNMPC(UnconstrParNMPC&&) noexcept = default;

  ///
  /// @brief Default move assign operator. 
  ///
  UnconstrParNMPC& operator=(UnconstrParNMPC&&) noexcept = default;

  ///
  /// @return const reference to the Robot model. 
  ///
  const Robot& robot() const {
    return robot_;
  }

  ///
  /// @return const reference to the cost function. 
  ///
  const std::shared_ptr<CostFunction>& cost() const {
    return cost_;
  }

  ///
  /// @return const reference to the constraints. 
  ///
  const std::shared_ptr<Constraints>& constraints() const {
    return constraints_;
  }

  ///
  /// @return Length of the horizon. 
  ///
  double T() const {
    return T_;
  }

  ///
  /// @return Number of the discretization grids on the horizon.
  ///
  int N() const {
    return N_;
  }

  ///
  /// @brief Overload operator[] to access the data as std::vector. 
  ///
  SplitUnconstrParNMPC& operator[] (const int i) {
    assert(i >= 0);
    assert(i < data.size());
    return data[i];
  }

  ///
  /// @brief const version of hybrid_container::operator[]. 
  ///
  const SplitUnconstrParNMPC& operator[] (const int i) const {
    assert(i >= 0);
    assert(i < data.size());
    return data[i];
  }

  ///
  /// @brief Split optimal control problem data for the time stages.
  ///
  std::vector<SplitUnconstrParNMPC> data;

  ///
  /// @brief Split optimal control problem data for the terminal stage.
  ///
  TerminalUnconstrParNMPC terminal;

  ///
  /// @brief Displays the optimal control problem onto a ostream.
  ///
  void disp(std::ostream& os) const {
    os << "UnconstrParNMPC: " << std::endl;
    os << "T: " << T_ << std::endl;
    os << "N: " << N_ << std::endl;
    os << robot_ << std::endl;
  }

  friend std::ostream& operator<<(std::ostream& os, 
                                  const UnconstrParNMPC& parnmpc) {
    parnmpc.disp(os);
    return os;
  }

private:
  Robot robot_;
  std::shared_ptr<CostFunction> cost_;
  std::shared_ptr<Constraints> constraints_;
  double T_;
  int N_;
};

} // namespace robotoc

#endif // ROBOTOC_UNCONSTR_PARNMPC_HPP_ 