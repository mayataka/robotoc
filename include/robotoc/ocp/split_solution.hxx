#ifndef ROBOTOC_SPLIT_SOLUTION_HXX_
#define ROBOTOC_SPLIT_SOLUTION_HXX_

#include "robotoc/ocp/split_solution.hpp"

#include <random>
#include <cassert>

namespace robotoc {

inline SplitSolution::SplitSolution(const Robot& robot) 
  : q(Eigen::VectorXd::Zero(robot.dimq())),
    v(Eigen::VectorXd::Zero(robot.dimv())),
    a(Eigen::VectorXd::Zero(robot.dimv())),
    u(Eigen::VectorXd::Zero(robot.dimu())),
    f(robot.maxPointContacts(), Eigen::Vector3d::Zero()),
    lmd(Eigen::VectorXd::Zero(robot.dimv())),
    gmm(Eigen::VectorXd::Zero(robot.dimv())),
    beta(Eigen::VectorXd::Zero(robot.dimv())),
    mu(robot.maxPointContacts(), Eigen::Vector3d::Zero()),
    nu_passive(Eigen::VectorXd::Zero(robot.dim_passive())),
    f_stack_(Eigen::VectorXd::Zero(robot.max_dimf())),
    mu_stack_(Eigen::VectorXd::Zero(robot.max_dimf())),
    xi_stack_(Eigen::VectorXd::Zero(robot.max_dimf())),
    has_floating_base_(robot.hasFloatingBase()),
    has_active_contacts_(false),
    has_active_impulse_(false),
    is_contact_active_(robot.maxPointContacts(), false),
    dimf_(0),
    dimi_(0) {
  if (robot.hasFloatingBase()) {
    q.coeffRef(6) = 1.0;
  }
}


inline SplitSolution::SplitSolution() 
  : q(),
    v(),
    a(),
    u(),
    f(),
    lmd(),
    gmm(),
    beta(),
    mu(),
    nu_passive(),
    f_stack_(),
    mu_stack_(),
    xi_stack_(),
    has_floating_base_(false),
    has_active_contacts_(false),
    has_active_impulse_(false),
    is_contact_active_(),
    dimf_(0),
    dimi_(0) {
}


inline SplitSolution::~SplitSolution() {
}


inline void SplitSolution::setContactStatus(
    const ContactStatus& contact_status) {
  assert(contact_status.maxPointContacts() == is_contact_active_.size());
  has_active_contacts_ = contact_status.hasActiveContacts();
  is_contact_active_ = contact_status.isContactActive();
  dimf_ = contact_status.dimf();
}


inline void SplitSolution::setContactStatus(const SplitSolution& other) {
  assert(other.isContactActive().size() == is_contact_active_.size());
  has_active_contacts_ = other.hasActiveContacts();
  is_contact_active_ = other.isContactActive();
  dimf_ = other.dimf();
}


inline void SplitSolution::setImpulseStatus(
    const ImpulseStatus& impulse_status) {
  has_active_impulse_ = impulse_status.hasActiveImpulse();
  dimi_ = impulse_status.dimf();
}


inline void SplitSolution::setImpulseStatus(const SplitSolution& other) {
  has_active_impulse_ = other.hasActiveImpulse();
  dimi_ = other.dimi();
}


inline void SplitSolution::setImpulseStatus() {
  has_active_impulse_ = false;
  dimi_ = 0;
}


inline Eigen::VectorBlock<Eigen::VectorXd> SplitSolution::f_stack() {
  return f_stack_.head(dimf_);
}


inline const Eigen::VectorBlock<const Eigen::VectorXd> 
SplitSolution::f_stack() const {
  return f_stack_.head(dimf_);
}


inline void SplitSolution::set_f_stack() {
  int contact_index = 0;
  int segment_start = 0;
  for (const auto is_contact_active : is_contact_active_) {
    if (is_contact_active) {
      f_stack_.template segment<3>(segment_start) = f[contact_index];
      segment_start += 3;
    }
    ++contact_index;
  }
}


inline void SplitSolution::set_f_vector() {
  int contact_index = 0;
  int segment_start = 0;
  for (const auto is_contact_active : is_contact_active_) {
    if (is_contact_active) {
      f[contact_index] = f_stack_.template segment<3>(segment_start);
      segment_start += 3;
    }
    ++contact_index;
  }
}


inline Eigen::VectorBlock<Eigen::VectorXd> SplitSolution::mu_stack() {
  return mu_stack_.head(dimf_);
}


inline const Eigen::VectorBlock<const Eigen::VectorXd> 
SplitSolution::mu_stack() const {
  return mu_stack_.head(dimf_);
}


inline void SplitSolution::set_mu_stack() {
  int contact_index = 0;
  int segment_start = 0;
  for (const auto is_contact_active : is_contact_active_) {
    if (is_contact_active) {
      mu_stack_.template segment<3>(segment_start) = mu[contact_index];
      segment_start += 3;
    }
    ++contact_index;
  }
}


inline void SplitSolution::set_mu_vector() {
  int contact_index = 0;
  int segment_start = 0;
  for (const auto is_contact_active : is_contact_active_) {
    if (is_contact_active) {
      mu[contact_index] = mu_stack_.template segment<3>(segment_start);
      segment_start += 3;
    }
    ++contact_index;
  }
}


inline Eigen::VectorBlock<Eigen::VectorXd> SplitSolution::xi_stack() {
  return xi_stack_.head(dimi_);
}


inline const Eigen::VectorBlock<const Eigen::VectorXd> 
SplitSolution::xi_stack() const {
  return xi_stack_.head(dimi_);
}


inline int SplitSolution::dimf() const {
  return dimf_;
}


inline int SplitSolution::dimi() const {
  return dimi_;
}


inline bool SplitSolution::hasActiveContacts() const {
  return has_active_contacts_;
}


inline bool SplitSolution::hasActiveImpulse() const {
  return has_active_impulse_;
}


inline bool SplitSolution::isContactActive(const int contact_index) const {
  assert(!is_contact_active_.empty());
  assert(contact_index >= 0);
  assert(contact_index < is_contact_active_.size());
  return is_contact_active_[contact_index];
}


inline std::vector<bool> SplitSolution::isContactActive() const {
  return is_contact_active_;
}


inline void SplitSolution::integrate(const Robot& robot, const double step_size, 
                                     const SplitDirection& d) {
  assert(f_stack().size() == d.df().size());
  assert(mu_stack().size() == d.dmu().size());
  robot.integrateConfiguration(d.dq(), step_size, q);
  v.noalias() += step_size * d.dv();
  a.noalias() += step_size * d.da();
  u.noalias() += step_size * d.du;
  lmd.noalias() += step_size * d.dlmd();
  gmm.noalias() += step_size * d.dgmm();
  beta.noalias() += step_size * d.dbeta();
  if (has_floating_base_) {
    nu_passive.noalias() += step_size * d.dnu_passive;
  }
  if (has_active_contacts_) {
    f_stack().noalias() += step_size * d.df();
    set_f_vector();
    mu_stack().noalias() += step_size * d.dmu();
    set_mu_vector();
  }
  if (has_active_impulse_) {
    assert(xi_stack().size() == d.dxi().size());
    xi_stack().noalias() += step_size * d.dxi();
  }
}


inline void SplitSolution::copyPrimal(const SplitSolution& other) {
  setContactStatus(other);
  q = other.q;
  v = other.v;
  a = other.a;
  u = other.u;
  for (int i=0; i<f.size(); ++i) {
    f[i] = other.f[i];
  }
  set_f_stack();
}


inline void SplitSolution::copyDual(const SplitSolution& other) {
  setContactStatus(other);
  setImpulseStatus(other);
  lmd = other.lmd;
  gmm = other.gmm;
  beta = other.beta;
  if (has_floating_base_) {
    nu_passive = other.nu_passive;
  }
  for (int i=0; i<f.size(); ++i) {
    mu[i] = other.mu[i];
  }
  set_mu_stack();
  if (has_active_impulse_) {
    xi_stack() = other.xi_stack();
  }
}

inline double SplitSolution::lagrangeMultiplierLinfNorm() const {
  const double lmd_linf = lmd.template lpNorm<Eigen::Infinity>();
  const double gmm_linf = gmm.template lpNorm<Eigen::Infinity>();
  const double beta_linf = beta.template lpNorm<Eigen::Infinity>();
  const double nu_passive_linf = (has_floating_base_ ? nu_passive.template lpNorm<Eigen::Infinity>() : 0);
  const double mu_linf = ((dimf_ > 0) ? mu_stack().template lpNorm<Eigen::Infinity>() : 0);
  const double xi_linf = ((dimi_ > 0) ? xi_stack().template lpNorm<Eigen::Infinity>() : 0);
  return std::max({lmd_linf, gmm_linf, beta_linf, nu_passive_linf, mu_linf, xi_linf});
}

inline bool SplitSolution::isApprox(const SplitSolution& other) const {
  if (!q.isApprox(other.q)) {
    return false;
  }
  if (!v.isApprox(other.v)) {
    return false;
  }
  if (!a.isApprox(other.a)) {
    return false;
  }
  if (!u.isApprox(other.u)) {
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
  if (has_floating_base_) {
    if (!nu_passive.isApprox(other.nu_passive)) {
      return false;
    }
  }
  if (has_active_contacts_) {
    if (!f_stack().isApprox(other.f_stack())) {
      return false;
    }
    if (!mu_stack().isApprox(other.mu_stack())) {
      return false;
    }
    for (int i=0; i<is_contact_active_.size(); ++i) {
      if (is_contact_active_[i]) {
        if (!other.isContactActive(i)) {
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
        if (other.isContactActive(i)) {
          return false;
        }
      }
    }
  }
  if (has_active_impulse_) {
    if (!xi_stack().isApprox(other.xi_stack())) {
      return false;
    }
  }
  return true;
}


inline void SplitSolution::setRandom(const Robot& robot) {
  q.setRandom(); 
  robot.normalizeConfiguration(q);
  v.setRandom();
  a.setRandom(); 
  u.setRandom(); 
  lmd.setRandom();
  gmm.setRandom(); 
  beta.setRandom(); 
  if (robot.hasFloatingBase()) {
    nu_passive.setRandom();
  }
}


inline void SplitSolution::setRandom(const Robot& robot, 
                                     const ContactStatus& contact_status) {
  setContactStatus(contact_status);
  setRandom(robot);
  if (contact_status.hasActiveContacts()) {
    f_stack().setRandom();
    mu_stack().setRandom();
    set_f_vector();
    set_mu_vector();
  }
}


inline void SplitSolution::setRandom(const Robot& robot, 
                                     const ImpulseStatus& impulse_status) {
  setImpulseStatus(impulse_status);
  setRandom(robot);
  if (impulse_status.hasActiveImpulse()) {
    xi_stack().setRandom();
  }
}


inline void SplitSolution::setRandom(const Robot& robot, 
                                     const ContactStatus& contact_status,
                                     const ImpulseStatus& impulse_status) {
  setRandom(robot, contact_status);
  setImpulseStatus(impulse_status);
  if (impulse_status.hasActiveImpulse()) {
    xi_stack().setRandom();
  }
}


inline SplitSolution SplitSolution::Random(const Robot& robot) {
  SplitSolution s(robot);
  s.setRandom(robot);
  return s;
}


inline SplitSolution SplitSolution::Random(
    const Robot& robot, const ContactStatus& contact_status) {
  SplitSolution s(robot);
  s.setRandom(robot, contact_status);
  return s;
}


inline SplitSolution SplitSolution::Random(
    const Robot& robot, const ImpulseStatus& impulse_status) {
  SplitSolution s(robot);
  s.setRandom(robot, impulse_status);
  return s;
}


inline SplitSolution SplitSolution::Random(
    const Robot& robot, const ContactStatus& contact_status, 
    const ImpulseStatus& impulse_status) {
  SplitSolution s(robot);
  s.setRandom(robot, contact_status, impulse_status);
  return s;
}

} // namespace robotoc 

#endif // ROBOTOC_SPLIT_SOLUTION_HXX_