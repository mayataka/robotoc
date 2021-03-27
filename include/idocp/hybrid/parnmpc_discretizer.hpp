#ifndef IDOCP_PARNMPC_DISCRETIZER_HPP_
#define IDOCP_PARNMPC_DISCRETIZER_HPP_

#include "idocp/hybrid/contact_sequence.hpp"

#include <vector>
#include <limits>

namespace idocp {


///
/// @class ParNMPCDiscretizer
/// @brief Discretizer of the hybrid optimal control problems for ParNMPC 
/// (backward Euler).
///
class ParNMPCDiscretizer {
public:
  ///
  /// @brief Constructor. 
  /// @param[in] T Length of the horizon.
  /// @param[in] N Ideal number of the discretization grids of the horizon. 
  /// Note that the actual number of the grids can differ from this value 
  /// depending on the discrete events.
  /// @param[in] max_num_events Maximum number of each discrete events 
  /// (impulse and lift). 
  ///
  ParNMPCDiscretizer(const double T, const int N, const int max_events);

  ///
  /// @brief Default constructor. 
  ///
  ParNMPCDiscretizer();

  ///
  /// @brief Destructor. 
  ///
  ~ParNMPCDiscretizer();

  ///
  /// @brief Default copy constructor. 
  ///
  ParNMPCDiscretizer(const ParNMPCDiscretizer&) = default;

  ///
  /// @brief Default copy assign operator. 
  ///
  ParNMPCDiscretizer& operator=(const ParNMPCDiscretizer&) = default;

  ///
  /// @brief Default move constructor. 
  ///
  ParNMPCDiscretizer(ParNMPCDiscretizer&&) noexcept = default;

  ///
  /// @brief Default move assign operator. 
  ///
  ParNMPCDiscretizer& operator=(ParNMPCDiscretizer&&) noexcept = default;

  ///
  /// @brief Discretizes the hybrid optimal control problem. 
  /// @param[in] contact_sequence Contact sequence.
  /// @param[in] t Initial time of the horizon.
  ///
  void discretizeOCP(const ContactSequence& contact_sequence, const double t);

  ///
  /// @return Number of the time stages on the horizon. 
  ///
  int N() const;

  ///
  /// @return Number of the impulse stages on the horizon. 
  ///
  int N_impulse() const;

  ///
  /// @return Number of the lift stages on the horizon. 
  ///
  int N_lift() const; 

  ///
  /// @return Total number of the time stages, impulse stages, auxiliary stages, 
  /// lift stages, and terminal stage on the horizon. 
  ///
  int N_all() const;

  ///
  /// @return Ideal number of the discretization grids on the horizon. 
  ///
  int N_ideal() const;

  ///
  /// @brief Returns the contact phase of the time stage. 
  /// @param[in] time_stage Time stage of interest. 
  /// @return Contact phase of the time stage. 
  ///
  int contactPhase(const int time_stage) const;

  ///
  /// @brief Returns the contact phase before the impulse. 
  /// @param[in] impulse_index Index of the impulse of interest. 
  /// @return Contact phase before the impulse. 
  ///
  int contactPhaseBeforeImpulse(const int impulse_index) const;

  ///
  /// @brief Returns the contact phase before the lift. 
  /// @param[in] impulse_index Index of the lift of interest. 
  /// @return Contact phase before the lift. 
  ///
  int contactPhaseBeforeLift(const int lift_index) const;

  ///
  /// @brief Returns the impulse index before the time stage. 
  /// @param[in] time_stage Time stage of interest. 
  /// @return Impulse index before the time stage. 
  ///
  int impulseIndexBeforeTimeStage(const int time_stage) const;

  ///
  /// @brief Returns the impulse index after the time stage. 
  /// @param[in] time_stage Time stage of interest. 
  /// @return Impulse index after the time stage. 
  ///
  int impulseIndexAfterTimeStage(const int time_stage) const;

  ///
  /// @brief Returns the lift index before the time stage. 
  /// @param[in] time_stage Time stage of interest. 
  /// @return Lift index before the time stage. 
  ///
  int liftIndexBeforeTimeStage(const int time_stage) const;

  ///
  /// @brief Returns the lift index after the time stage. 
  /// @param[in] time_stage Time stage of interest. 
  /// @return Lift index after the time stage. 
  ///
  int liftIndexAfterTimeStage(const int time_stage) const;

  ///
  /// @brief Returns the time stage before the impulse. 
  /// @param[in] impulse_index Index of the impulse of interest. 
  /// @return Time stage before the impulse. 
  ///
  int timeStageBeforeImpulse(const int impulse_index) const;

  ///
  /// @brief Returns the time stage after the impulse. 
  /// @param[in] impulse_index Index of the impulse of interest. 
  /// @return Time stage after the impulse. 
  ///
  int timeStageAfterImpulse(const int impulse_index) const;

  ///
  /// @brief Returns the time stage before the lift. 
  /// @param[in] impulse_index Index of the lift of interest. 
  /// @return Time stage before the lift. 
  ///
  int timeStageBeforeLift(const int lift_index) const;

  ///
  /// @brief Returns the time stage after the lift. 
  /// @param[in] impulse_index Index of the lift of interest. 
  /// @return Time stage after the lift. 
  ///
  int timeStageAfterLift(const int lift_index) const;

  ///
  /// @brief Checks wheather the time stage is just before the impulse or not. 
  /// @param[in] time_stage Time stage of interest. 
  /// @return true if the time stage is just before the impulse. false if not.
  ///
  bool isTimeStageBeforeImpulse(const int time_stage) const;

  ///
  /// @brief Checks wheather the time stage is just after the impulse or not. 
  /// @param[in] time_stage Time stage of interest. 
  /// @return true if the time stage is just after the impulse. false if not.
  ///
  bool isTimeStageAfterImpulse(const int time_stage) const;

  ///
  /// @brief Checks wheather the time stage is just before the lift or not. 
  /// @param[in] time_stage Time stage of interest. 
  /// @return true if the time stage is just before the lift. false if not.
  ///
  bool isTimeStageBeforeLift(const int time_stage) const;

  ///
  /// @brief Checks wheather the time stage is just after the lift or not. 
  /// @param[in] time_stage Time stage of interest. 
  /// @return true if the time stage is just after the lift. false if not.
  ///
  bool isTimeStageAfterLift(const int time_stage) const;

  ///
  /// @brief Returns the time of the time stage. 
  /// @param[in] time_stage Time stage of interest. 
  /// @return Time of the time stage of interest.
  ///
  double t(const int time_stage) const;

  ///
  /// @brief Returns the time of the impulse. 
  /// @param[in] impulse_index Index of impulse of interest. 
  /// @return Time of the impulse of interest.
  ///
  double t_impulse(const int impulse_index) const;

  ///
  /// @brief Returns the time of the lift. 
  /// @param[in] impulse_index Index of lift of interest. 
  /// @return Time of the lift of interest.
  ///
  double t_lift(const int lift_index) const;

  ///
  /// @brief Returns the time step of the time stage. 
  /// @param[in] time_stage Time stage of interest. 
  /// @return Time step of the time stage of interest.
  ///
  double dt(const int time_stage) const;

  ///
  /// @brief Returns the time step of the auxiliary stage of the impulse. 
  /// @param[in] impulse_index Index of impulse of interest. 
  /// @return Time step of the auxiliary stage of the impulse.
  ///
  double dt_aux(const int impulse_index) const;

  ///
  /// @brief Returns the time step of the auxiliary stage of the lift. 
  /// @param[in] lift_index Index of lift of interest. 
  /// @return Time step of the auxiliary stage of the lift.
  ///
  double dt_lift(const int lift_index) const;

  ///
  /// @brief Checks wheather the optimal control problem is well-defined. 
  /// @return true if the optimal control problem is well-defined. false if not.
  ///
  bool isWellDefined() const;

  ///
  /// @brief Shows the information of the discretized optimal control problem 
  /// into console. 
  ///
  void showInfo() const;

private:
  double T_, dt_ideal_, max_dt_;
  int N_, N_ideal_, N_impulse_, N_lift_, max_events_;
  std::vector<int> contact_phase_index_from_time_stage_, 
                   impulse_index_before_time_stage_, 
                   lift_index_before_time_stage_, time_stage_after_impulse_, 
                   time_stage_after_lift_;
  std::vector<bool> is_time_stage_after_impulse_, is_time_stage_after_lift_;
  std::vector<double> t_, t_impulse_, t_lift_, dt_, dt_aux_, dt_lift_;

  static constexpr double min_dt_ 
      = std::sqrt(std::numeric_limits<double>::epsilon());

  void countDiscreteEvents(const ContactSequence& contact_sequence, 
                           const double t);

  void countTimeSteps(const double t);

  void countTimeStages();

  void countContactPhase();

};

} // namespace idocp

#include "idocp/hybrid/parnmpc_discretizer.hxx"

#endif // IDOCP_PARNMPC_DISCRETIZER_HPP_ 