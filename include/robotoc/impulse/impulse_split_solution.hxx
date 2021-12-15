#ifndef ROBOTOC_IMPULSE_SPLIT_SOLUTION_HXX_
#define ROBOTOC_IMPULSE_SPLIT_SOLUTION_HXX_

#include "robotoc/impulse/impulse_split_solution.hpp"

#include <random>

namespace robotoc {

inline ImpulseSplitSolution::ImpulseSplitSolution(const Robot& robot) 
  : q(Eigen::VectorXd::Zero(robot.dimq())),
    v(Eigen::VectorXd::Zero(robot.dimv())),
    dv(Eigen::VectorXd::Zero(robot.dimv())),
    f(robot.maxPointContacts(), Eigen::Vector3d::Zero()),
    lmd(Eigen::VectorXd::Zero(robot.dimv())),
    gmm(Eigen::VectorXd::Zero(robot.dimv())),
    beta(Eigen::VectorXd::Zero(robot.dimv())),
    mu(robot.maxPointContacts(), Eigen::Vector3d::Zero()),
    f_stack_(Eigen::VectorXd::Zero(robot.max_dimf())),
    mu_stack_(Eigen::VectorXd::Zero(robot.max_dimf())),
    is_impulse_active_(robot.maxPointContacts(), false),
    dimi_(0) {
  if (robot.hasFloatingBase()) {
    q.coeffRef(6) = 1.0;
  }
}


inline ImpulseSplitSolution::ImpulseSplitSolution() 
  : q(),
    v(),
    dv(),
    f(),
    lmd(),
    gmm(),
    beta(),
    mu(),
    f_stack_(),
    mu_stack_(),
    is_impulse_active_(),
    dimi_(0) {
}


inline ImpulseSplitSolution::~ImpulseSplitSolution() {
}


inline void ImpulseSplitSolution::setImpulseStatus(
    const ImpulseStatus& impulse_status) {
  assert(impulse_status.maxPointContacts() == is_impulse_active_.size());
  is_impulse_active_ = impulse_status.isImpulseActive();
  dimi_ = impulse_status.dimi();
}


inline void ImpulseSplitSolution::setImpulseStatus(
    const ImpulseSplitSolution& other) {
  assert(other.isImpulseActive().size() == is_impulse_active_.size());
  is_impulse_active_ = other.isImpulseActive();
  dimi_ = other.dimi();
}


inline Eigen::VectorBlock<Eigen::VectorXd> ImpulseSplitSolution::f_stack() {
  return f_stack_.head(dimi_);
}


inline const Eigen::VectorBlock<const Eigen::VectorXd> 
ImpulseSplitSolution::f_stack() const {
  return f_stack_.head(dimi_);
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
  return mu_stack_.head(dimi_);
}


inline const Eigen::VectorBlock<const Eigen::VectorXd> 
ImpulseSplitSolution::mu_stack() const {
  return mu_stack_.head(dimi_);
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


inline int ImpulseSplitSolution::dimi() const {
  return dimi_;
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


inline void ImpulseSplitSolution::integrate(const Robot& robot, 
                                            const double step_size, 
                                            const ImpulseSplitDirection& d) {
  assert(f_stack().size() == d.df().size());
  assert(mu_stack().size() == d.dmu().size());
  robot.integrateConfiguration(d.dq(), step_size, q);
  v.noalias() += step_size * d.dv();
  dv.noalias() += step_size * d.ddv();
  f_stack().noalias() += step_size * d.df();
  set_f_vector();
  lmd.noalias() += step_size * d.dlmd();
  gmm.noalias() += step_size * d.dgmm();
  beta.noalias() += step_size * d.dbeta();
  mu_stack().noalias() += step_size * d.dmu();
  set_mu_vector();
}


inline void ImpulseSplitSolution::copyPrimal(const ImpulseSplitSolution& other) {
  setImpulseStatus(other);
  q = other.q;
  v = other.v;
  dv = other.dv;
  for (int i=0; i<f.size(); ++i) {
    f[i] = other.f[i];
  }
  set_f_stack();
}


inline void ImpulseSplitSolution::copyDual(const ImpulseSplitSolution& other) {
  setImpulseStatus(other);
  lmd = other.lmd;
  gmm = other.gmm;
  beta = other.beta;
  for (int i=0; i<f.size(); ++i) {
    mu[i] = other.mu[i];
  }
  set_mu_stack();
}

inline double ImpulseSplitSolution::lagrangeMultiplierLinfNorm() const {
  const double lmd_linf = lmd.template lpNorm<Eigen::Infinity>();
  const double gmm_linf = gmm.template lpNorm<Eigen::Infinity>();
  const double beta_linf = beta.template lpNorm<Eigen::Infinity>();
  const double mu_linf = ((dimi_ > 0) ? mu_stack().template lpNorm<Eigen::Infinity>() : 0);
  return std::max({lmd_linf, gmm_linf, beta_linf, mu_linf});
}

inline bool ImpulseSplitSolution::isApprox(
    const ImpulseSplitSolution& other) const {
  if (!q.isApprox(other.q)) {
    return false;
  }
  if (!v.isApprox(other.v)) {
    return false;
  }
  if (!dv.isApprox(other.dv)) {
    return false;
  }
  if (!f_stack().isApprox(other.f_stack())) {
    return false;
  }
  if (!lmd.isApprox(other.lmd)) {
    return false;
  }
  if (!gmm.isApprox(other.gmm)) {
    return false;
  }
  if (!beta.isApprox(other.beta)) {
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
  q.setRandom(); 
  robot.normalizeConfiguration(q);
  v.setRandom();
  dv.setRandom(); 
  lmd.setRandom();
  gmm.setRandom(); 
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

} // namespace robotoc 

#endif // ROBOTOC_IMPULSE_SPLIT_SOLUTION_HXX_ 