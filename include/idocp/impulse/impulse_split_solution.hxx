#ifndef IDOCP_IMPULSE_SPLIT_SOLUTION_HXX_
#define IDOCP_IMPULSE_SPLIT_SOLUTION_HXX_

#include "idocp/impulse/impulse_split_solution.hpp"

namespace idocp {

inline ImpulseSplitSolution::ImpulseSplitSolution(const Robot& robot) 
  : lmd(Eigen::VectorXd::Zero(robot.dimv())),
    gmm(Eigen::VectorXd::Zero(robot.dimv())),
    mu_contact_velocity(robot.max_point_contacts(), Eigen::Vector3d::Zero()),
    mu_contact_position(robot.max_point_contacts(), Eigen::Vector3d::Zero()),
    dv(Eigen::VectorXd::Zero(robot.dimv())),
    f(robot.max_point_contacts(), Eigen::Vector3d::Zero()),
    q(Eigen::VectorXd::Zero(robot.dimq())),
    v(Eigen::VectorXd::Zero(robot.dimv())),
    beta(Eigen::VectorXd::Zero(robot.dimv())),
    mu_stack_(Eigen::VectorXd::Zero(2*robot.max_dimf())),
    f_stack_(Eigen::VectorXd::Zero(robot.max_dimf())),
    is_contact_active_(robot.max_point_contacts(), false),
    dimf_(0),
    dimc_(0) {
  robot.normalizeConfiguration(q);
}


inline ImpulseSplitSolution::ImpulseSplitSolution() 
  : lmd(),
    gmm(),
    mu_contact_velocity(),
    mu_contact_position(),
    dv(),
    f(),
    q(),
    v(),
    beta(),
    mu_stack_(),
    f_stack_(),
    is_contact_active_(),
    dimf_(0),
    dimc_(0) {
}


inline ImpulseSplitSolution::~ImpulseSplitSolution() {
}


inline void ImpulseSplitSolution::setContactStatus(
    const ContactStatus& contact_status) {
  is_contact_active_ = contact_status.isContactActive();
  dimf_ = contact_status.dimf();
}


inline Eigen::VectorBlock<Eigen::VectorXd> ImpulseSplitSolution::mu_stack() {
  return mu_stack_.head(dimc_);
}


inline const Eigen::VectorBlock<const Eigen::VectorXd> 
ImpulseSplitSolution::mu_stack() const {
  return mu_stack_.head(dimc_);
}


inline Eigen::VectorBlock<Eigen::VectorXd> ImpulseSplitSolution::f_stack() {
  return f_stack_.head(dimf_);
}


inline const Eigen::VectorBlock<const Eigen::VectorXd> 
ImpulseSplitSolution::f_stack() const {
  return f_stack_.head(dimf_);
}


inline void ImpulseSplitSolution::set_mu_stack() {
  int contact_index = 0;
  int segment_start = 0;
  for (const auto is_contact_active : is_contact_active_) {
    if (is_contact_active) {
      mu_stack_.template segment<3>(segment_start) 
          = mu_contact_velocity[contact_index];
      segment_start += 3;
    }
    ++contact_index;
  }
  contact_index = 0;
  for (const auto is_contact_active : is_contact_active_) {
    if (is_contact_active) {
      mu_stack_.template segment<3>(segment_start) 
        = mu_contact_position[contact_index];
      segment_start += 3;
    }
    ++contact_index;
  }
}


inline void ImpulseSplitSolution::set_mu_contact() {
  int contact_index = 0;
  int segment_start = 0;
  for (const auto is_contact_active : is_contact_active_) {
    if (is_contact_active) {
      mu_contact_velocity[contact_index]
        = mu_stack_.template segment<3>(segment_start);
      segment_start += 3;
    }
    ++contact_index;
  }
  contact_index = 0;
  for (const auto is_contact_active : is_contact_active_) {
    if (is_contact_active) {
      mu_contact_position[contact_index]
        = mu_stack_.template segment<3>(segment_start);
      segment_start += 3;
    }
    ++contact_index;
  }
}


inline void ImpulseSplitSolution::set_f_stack() {
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


inline void ImpulseSplitSolution::set_f() {
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


inline bool ImpulseSplitSolution::isContactActive(
    const int contact_index) const {
  return is_contact_active_[contact_index];
}


inline int ImpulseSplitSolution::num_active_contacts() const {
  return is_contact_active_.size();
}


inline int ImpulseSplitSolution::dimc() const {
  return dimc_;
}


inline int ImpulseSplitSolution::dimf() const {
  return dimf_;
}


inline ImpulseSplitSolution ImpulseSplitSolution::Random(
    const Robot& robot, const ContactStatus& contact_status) {
  ImpulseSplitSolution s(robot);
  s.setContactStatus(contact_status);
  s.lmd = Eigen::VectorXd::Random(robot.dimv());
  s.gmm = Eigen::VectorXd::Random(robot.dimv());
  s.mu_stack() = Eigen::VectorXd::Random(s.dimc());
  s.set_mu_contact();
  s.dv = Eigen::VectorXd::Random(robot.dimv());
  s.f_stack() = Eigen::VectorXd::Random(s.dimf());
  s.set_f();
  s.q = Eigen::VectorXd::Random(robot.dimq());
  robot.normalizeConfiguration(s.q);
  s.v = Eigen::VectorXd::Random(robot.dimv());
  s.beta = Eigen::VectorXd::Random(robot.dimv());
  return s;
}

} // namespace idocp 

#endif // IDOCP_IMPULSE_SPLIT_SOLUTION_HXX_ 