#ifndef IDOCP_SPLIT_TEMPORARY_SOLUTION_HXX_
#define IDOCP_SPLIT_TEMPORARY_SOLUTION_HXX_

#include "idocp/line_search/split_temporary_solution.hpp"

#include <assert.h>

namespace idocp {

inline SplitTemporarySolution::SplitTemporarySolution(const Robot& robot) 
  : s_tmp_(SplitSolution(robot)),
    q_next_tmp_(Eigen::VectorXd::Zero(robot.dimq())),
    v_next_tmp_(Eigen::VectorXd::Zero(robot.dimv())) {
}


inline SplitTemporarySolution::SplitTemporarySolution() 
  : s_tmp_(),
    q_next_tmp_(),
    v_next_tmp_() {
}


inline SplitTemporarySolution::~SplitTemporarySolution() {
}


inline void SplitTemporarySolution::setTemporarySolution(
    const Robot& robot, const ContactStatus& contact_status, 
    const double step_size, const SplitSolution& s, const SplitDirection& d, 
    const SplitSolution& s_next, const SplitDirection& d_next) {
  s_tmp_.setContactStatus(contact_status);
  robot.integrateConfiguration(s.q, d.dq(), step_size, s_tmp_.q);
  s_tmp_.v = s.v + step_size * d.dv();
  s_tmp_.a = s.a + step_size * d.da();
  if (contact_status.hasActiveContacts()) {
    s_tmp_.f_stack() = s.f_stack() + step_size * d.df();
  }
  s_tmp_.u = s.u + step_size * d.du();
  if (robot.has_floating_base()) {
    s_tmp_.u_passive = s.u_passive + step_size * d.du_passive;
  }
  robot.integrateConfiguration(s_next.q, d_next.dq(), step_size, q_next_tmp_);
  v_next_tmp_ = s_next.v + step_size * d_next.dv();
}


inline void SplitTemporarySolution::setTemporarySolution(
    const Robot& robot, const double step_size, const SplitSolution& s, 
    const SplitDirection& d) {
  robot.integrateConfiguration(s.q, d.dq(), step_size, s_tmp_.q);
  s_tmp_.v = s.v + step_size * d.dv();
}


inline const SplitSolution& SplitTemporarySolution::splitSolution() const {
  return s_tmp_;
}


inline const Eigen::VectorXd& SplitTemporarySolution::q_next() const {
  return q_next_tmp_;
}


inline const Eigen::VectorXd& SplitTemporarySolution::v_next() const {
  return v_next_tmp_;
}


} // namespace idocp


#endif // IDOCP_SPLIT_TEMPORARY_SOLUTION_HXX_