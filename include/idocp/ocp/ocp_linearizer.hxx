#ifndef IDOCP_OCP_LINEARIZER_HXX_ 
#define IDOCP_OCP_LINEARIZER_HXX_

#include "idocp/ocp/ocp_linearizer.hpp"

#include <omp.h>
#include <stdexcept>
#include <cassert>

namespace idocp {

inline OCPLinearizer::OCPLinearizer(
    const Robot& robot, const std::shared_ptr<CostFunction>& cost, 
    const std::shared_ptr<Constraints>& constraints, const double T, 
    const int N, const int max_num_impulse, const int num_proc) 
  : split_ocps_(N, SplitOCP(robot, cost, constraints), 
                max_num_impulse, SplitImpulseOCP(robot, 
                                                 cost->getImpulseCostFunction(), 
                                                 constraints->getImpulseConstraints())),
    terminal_ocp_(robot, cost, constraints),
    T_(T),
    dtau_(T/N),
    N_(N),
    num_proc_(num_proc) {
  try {
    if (T <= 0) {
      throw std::out_of_range("invalid value: T must be positive!");
    }
    if (N <= 0) {
      throw std::out_of_range("invalid value: N must be positive!");
    }
    if (max_num_impulse < 0) {
      throw std::out_of_range("invalid value: max_num_impulse must be non-negative!");
    }
    if (num_proc <= 0) {
      throw std::out_of_range("invalid value: num_proc must be positive!");
    }
  }
  catch(const std::exception& e) {
    std::cerr << e.what() << '\n';
    std::exit(EXIT_FAILURE);
  }
}


inline OCPLinearizer::OCPLinearizer()
  : split_ocps_(),
    terminal_ocp_(),
    T_(0),
    dtau_(0),
    N_(0),
    num_proc_(0) {
}


inline OCPLinearizer::~OCPLinearizer() {
}


inline void OCPLinearizer::linearizeInitialState(const Robot& robot,
                                                 const double t, 
                                                 const Eigen::VectorXd& q, 
                                                 const Eigen::VectorXd& v,
                                                 const HybridSolution& s, 
                                                 HybridDirection& d) const {
  robot.subtractConfiguration(q, s[0].q, d[0].dq());
  d[0].dv() = v - s[0].v;
}


inline void OCPLinearizer::updateSolution(const std::vector<Robot>& robots, 
                                          const HybridSolution& s, 
                                          const double primal_step_size, 
                                          const double dual_step_size, 
                                          HybridDirection& d) {
}


inline const Eigen::VectorXd& OCPLinearizer::q_prev(
    const ContactSequence& contact_sequence, const HybridSolution& s,
    const int time_stage) const {
  assert(time_stage >= 1);
  assert(time_stage <= N_);
  if (contact_sequence.existImpulseStage(time_stage-1)) {
    return s.aux[contact_sequence.impulseIndex(time_stage-1)].q;
  }
  else if (contact_sequence.existLiftStage(time_stage-1)) {
    return s.lift[contact_sequence.liftIndex(time_stage-1)].q;
  }
  else {
    return s[time_stage-1].q;
  }
}

} // namespace idocp 

#endif // IDOCP_OCP_LINEARIZER_HXX_ 