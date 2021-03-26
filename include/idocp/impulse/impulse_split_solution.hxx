#ifndef IDOCP_IMPULSE_SPLIT_SOLUTION_HXX_
#define IDOCP_IMPULSE_SPLIT_SOLUTION_HXX_

#include "idocp/impulse/impulse_split_solution.hpp"

#include <random>

namespace idocp {

inline ImpulseSplitSolution::ImpulseSplitSolution(const Robot& robot) 
  : lmd(Eigen::VectorXd::Zero(robot.dimv())),
    gmm(Eigen::VectorXd::Zero(robot.dimv())),
    q(Eigen::VectorXd::Zero(robot.dimq())),
    v(Eigen::VectorXd::Zero(robot.dimv())),
    dv(Eigen::VectorXd::Zero(robot.dimv())),
    f(robot.maxPointContacts(), Eigen::Vector3d::Zero()),
    beta(Eigen::VectorXd::Zero(robot.dimv())),
    mu(robot.maxPointContacts(), Eigen::Vector3d::Zero()),
    f_stack_(Eigen::VectorXd::Zero(robot.max_dimf())),
    mu_stack_(Eigen::VectorXd::Zero(robot.max_dimf())),
    has_floating_base_(robot.hasFloatingBase()),
    is_impulse_active_(robot.maxPointContacts(), false),
    dimf_(0) {
  robot.normalizeConfiguration(q);
}


inline ImpulseSplitSolution::ImpulseSplitSolution() 
  : lmd(),
    gmm(),
    q(),
    v(),
    dv(),
    f(),
    beta(),
    mu(),
    f_stack_(),
    mu_stack_(),
    has_floating_base_(false),
    is_impulse_active_(),
    dimf_(0) {
}


inline ImpulseSplitSolution::~ImpulseSplitSolution() {
}


inline void ImpulseSplitSolution::setImpulseStatus(
    const ImpulseStatus& impulse_status) {
  assert(impulse_status.maxPointContacts() == is_impulse_active_.size());
  is_impulse_active_ = impulse_status.isImpulseActive();
  dimf_ = impulse_status.dimf();
}


inline void ImpulseSplitSolution::setImpulseStatus(
    const ImpulseSplitSolution& other) {
  assert(other.isImpulseActive().size() == is_impulse_active_.size());
  is_impulse_active_ = other.isImpulseActive();
  dimf_ = other.dimf();
}


inline Eigen::VectorBlock<Eigen::VectorXd> ImpulseSplitSolution::f_stack() {
  return f_stack_.head(dimf_);
}


inline const Eigen::VectorBlock<const Eigen::VectorXd> 
ImpulseSplitSolution::f_stack() const {
  return f_stack_.head(dimf_);
}


inline void ImpulseSplitSolution::set_f_stack() {
  int contact_index = 0;
  int segment_start = 0;
  for (const auto is_impulse_active : is_impulse_active_) {
    if (is_impulse_active) {
      f_stack_.template segment<3>(segment_start) = f[contact_index];
      segment_start += 3;
    }
    ++contact_index;
  }
}


inline void ImpulseSplitSolution::set_f_vector() {
  int contact_index = 0;
  int segment_start = 0;
  for (const auto is_impulse_active : is_impulse_active_) {
    if (is_impulse_active) {
      f[contact_index] = f_stack_.template segment<3>(segment_start);
      segment_start += 3;
    }
    ++contact_index;
  }
}


inline Eigen::VectorBlock<Eigen::VectorXd> ImpulseSplitSolution::mu_stack() {
  return mu_stack_.head(dimf_);
}


inline const Eigen::VectorBlock<const Eigen::VectorXd> 
ImpulseSplitSolution::mu_stack() const {
  return mu_stack_.head(dimf_);
}


inline void ImpulseSplitSolution::set_mu_stack() {
  int contact_index = 0;
  int segment_start = 0;
  for (const auto is_impulse_active : is_impulse_active_) {
    if (is_impulse_active) {
      mu_stack_.template segment<3>(segment_start) = mu[contact_index];
      segment_start += 3;
    }
    ++contact_index;
  }
}


inline void ImpulseSplitSolution::set_mu_vector() {
  int contact_index = 0;
  int segment_start = 0;
  for (const auto is_impulse_active : is_impulse_active_) {
    if (is_impulse_active) {
      mu[contact_index] = mu_stack_.template segment<3>(segment_start);
      segment_start += 3;
    }
    ++contact_index;
  }
}


inline int ImpulseSplitSolution::dimf() const {
  return dimf_;
}


inline bool ImpulseSplitSolution::isImpulseActive(
    const int contact_index) const {
  assert(!is_impulse_active_.empty());
  assert(contact_index >= 0);
  assert(contact_index < is_impulse_active_.size());
  return is_impulse_active_[contact_index];
}


inline std::vector<bool> ImpulseSplitSolution::isImpulseActive() const {
  return is_impulse_active_;
}


inline void ImpulseSplitSolution::copy(const ImpulseSplitSolution& other) {
  setImpulseStatus(other);
  lmd        = other.lmd;
  gmm        = other.gmm;
  q          = other.q;
  v          = other.v;
  dv         = other.dv;
  beta       = other.beta;
  f_stack()  = other.f_stack();
  mu_stack() = other.mu_stack();
  set_f_vector();
  set_mu_vector();
}


inline void ImpulseSplitSolution::copyPartial(const SplitSolution& s) {
  lmd        = s.lmd;
  gmm        = s.gmm;
  q          = s.q;
  v          = s.v;
  dv         = s.a;
  beta       = s.beta;
  const int fsize = f.size();
  for (int i=0; i<fsize; ++i) {
    f[i] = s.f[i];
  }
  for (int i=0; i<fsize; ++i) {
    mu[i] = s.mu[i];
  }
}


inline void ImpulseSplitSolution::integrate(const Robot& robot, 
                                            const double step_size, 
                                            const ImpulseSplitDirection& d) {
  lmd.noalias() += step_size * d.dlmd();
  gmm.noalias() += step_size * d.dgmm();
  robot.integrateConfiguration(d.dq(), step_size, q);
  v.noalias() += step_size * d.dv();
  dv.noalias() += step_size * d.ddv();
  beta.noalias() += step_size * d.dbeta();
  assert(f_stack().size() == d.df().size());
  f_stack().noalias() += step_size * d.df();
  set_f_vector();
  assert(mu_stack().size() == d.dmu().size());
  mu_stack().noalias() += step_size * d.dmu();
  set_mu_vector();
}


inline bool ImpulseSplitSolution::isApprox(
    const ImpulseSplitSolution& other) const {
  if (!lmd.isApprox(other.lmd)) {
    return false;
  }
  if (!gmm.isApprox(other.gmm)) {
    return false;
  }
  if (!q.isApprox(other.q)) {
    return false;
  }
  if (!v.isApprox(other.v)) {
    return false;
  }
  if (!dv.isApprox(other.dv)) {
    return false;
  }
  if (!beta.isApprox(other.beta)) {
    return false;
  }
  if (!f_stack().isApprox(other.f_stack())) {
    return false;
  }
  if (!mu_stack().isApprox(other.mu_stack())) {
    return false;
  }
  for (int i=0; i<is_impulse_active_.size(); ++i) {
    if (is_impulse_active_[i]) {
      if (!other.isImpulseActive(i)) {
        return false;
      }
      if (!f[i].isApprox(other.f[i])) {
        return false;
      }
      if (!mu[i].isApprox(other.mu[i])) {
        return false;
      }
    }
    else {
      if (other.isImpulseActive(i)) {
        return false;
      }
    }
  }
  return true;
}


inline void ImpulseSplitSolution::setRandom(const Robot& robot) {
  lmd.setRandom();
  gmm.setRandom(); 
  q.setRandom(); 
  robot.normalizeConfiguration(q);
  v.setRandom();
  dv.setRandom(); 
  beta.setRandom(); 
}


inline void ImpulseSplitSolution::setRandom(
    const Robot& robot, const ImpulseStatus& impulse_status) {
  setImpulseStatus(impulse_status);
  setRandom(robot);
  if (impulse_status.hasActiveImpulse()) {
    f_stack().setRandom();
    mu_stack().setRandom();
    set_f_vector();
    set_mu_vector();
  }
}


inline ImpulseSplitSolution ImpulseSplitSolution::Random(const Robot& robot) {
  ImpulseSplitSolution s(robot);
  s.setRandom(robot);
  return s;
}


inline ImpulseSplitSolution ImpulseSplitSolution::Random(
    const Robot& robot, const ImpulseStatus& impulse_status) {
  ImpulseSplitSolution s(robot);
  s.setRandom(robot, impulse_status);
  return s;
}

} // namespace idocp 

#endif // IDOCP_IMPULSE_SPLIT_SOLUTION_HXX_ 