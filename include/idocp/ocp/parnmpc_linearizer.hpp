#ifndef IDOCP_PARNMPC_LINEARIZER_HPP_
#define IDOCP_PARNMPC_LINEARIZER_HPP_

#include <vector>
#include <memory>
#include <cassert>

#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/hybrid/hybrid_container.hpp"
#include "idocp/hybrid/contact_sequence.hpp"
#include "idocp/hybrid/parnmpc_discretizer.hpp"


namespace idocp {

///
/// @class ParNMPCLinearizer
/// @brief Linearize of the optimal control problem. 
///
class ParNMPCLinearizer {
public:
  ///
  /// @brief Construct optimal control problem solver.
  /// @param[in] N Number of discretization of the horizon. Must be more than 1. 
  /// @param[in] max_num_impulse Maximum number of the impulse on the horizon. 
  /// Must be non-negative. 
  /// @param[in] nthreads Number of the threads in solving the optimal control 
  /// problem. Must be positive. 
  ///
  ParNMPCLinearizer(const int N, const int max_num_impulse, const int nthreads);

  ///
  /// @brief Default constructor. 
  ///
  ParNMPCLinearizer();

  ///
  /// @brief Destructor. 
  ///
  ~ParNMPCLinearizer();

  ///
  /// @brief Default copy constructor. 
  ///
  ParNMPCLinearizer(const ParNMPCLinearizer&) = default;

  ///
  /// @brief Default copy assign operator. 
  ///
  ParNMPCLinearizer& operator=(const ParNMPCLinearizer&) = default;

  ///
  /// @brief Default move constructor. 
  ///
  ParNMPCLinearizer(ParNMPCLinearizer&&) noexcept = default;

  ///
  /// @brief Default move assign operator. 
  ///
  ParNMPCLinearizer& operator=(ParNMPCLinearizer&&) noexcept = default;

  void initConstraints(ParNMPC& parnmpc, std::vector<Robot>& robots,
                       const ContactSequence& contact_sequence, 
                       const Solution& s) const;

  void computeKKTResidual(ParNMPC& parnmpc, std::vector<Robot>& robots, 
                          const ContactSequence& contact_sequence,
                          const Eigen::VectorXd& q, const Eigen::VectorXd& v, 
                          const Solution& s, KKTMatrix& kkt_matrix, 
                          KKTResidual& kkt_residual) const;

  double KKTError(const ParNMPC& parnmpc, const KKTResidual& kkt_residual);

  void integrateSolution(ParNMPC& parnmpc, const std::vector<Robot>& robots,
                         const KKTMatrix& kkt_matrix,
                         const KKTResidual& kkt_residual,
                         const double primal_step_size,
                         const double dual_step_size,
                         Direction& d, Solution& s) const;

  static const Eigen::VectorXd& q_prev(const ParNMPCDiscretizer& discretizer, 
                                       const Eigen::VectorXd& q, 
                                       const Solution& s, 
                                       const int time_stage) {
    assert(time_stage >= 0);
    assert(time_stage < discretizer.N());
    if (discretizer.isTimeStageAfterImpulse(time_stage)) {
      return s.impulse[discretizer.impulseIndexBeforeTimeStage(time_stage)].q;
    }
    else if (discretizer.isTimeStageAfterLift(time_stage)) {
      return s.lift[discretizer.liftIndexBeforeTimeStage(time_stage)].q;
    }
    else if (time_stage > 0) {
      return s[time_stage-1].q;
    }
    else { 
      assert(time_stage == 0);
      return q;
    }
  }

  static const Eigen::VectorXd& v_prev(const ParNMPCDiscretizer& discretizer, 
                                       const Eigen::VectorXd& v, 
                                       const Solution& s, 
                                       const int time_stage) {
    assert(time_stage >= 0);
    assert(time_stage < discretizer.N());
    if (discretizer.isTimeStageAfterImpulse(time_stage)) {
      return s.impulse[discretizer.impulseIndexBeforeTimeStage(time_stage)].v;
    }
    else if (discretizer.isTimeStageAfterLift(time_stage)) {
      return s.lift[discretizer.liftIndexBeforeTimeStage(time_stage)].v;
    }
    else if (time_stage > 0) {
      return s[time_stage-1].v;
    }
    else { 
      assert(time_stage == 0);
      return v;
    }
  }

private:
  int max_num_impulse_, nthreads_;
  Eigen::VectorXd kkt_error_;

};

} // namespace idocp 

#endif // IDOCP_PARNMPC_LINEARIZER_HPP_ 