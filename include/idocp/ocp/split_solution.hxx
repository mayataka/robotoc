#ifndef IDOCP_SPLIT_SOLUTION_HXX_
#define IDOCP_SPLIT_SOLUTION_HXX_

#include "idocp/ocp/split_solution.hpp"

namespace idocp {

inline SplitSolution::SplitSolution(const Robot& robot) 
  : lmd(Eigen::VectorXd::Zero(robot.dimv())),
    gmm(Eigen::VectorXd::Zero(robot.dimv())),
    mu_contact(robot.max_point_contacts(), Eigen::Vector3d::Zero()),
    a(Eigen::VectorXd::Zero(robot.dimv())),
    f(robot.max_point_contacts(), Eigen::Vector3d::Zero()),
    q(Eigen::VectorXd::Zero(robot.dimq())),
    v(Eigen::VectorXd::Zero(robot.dimv())),
    u(Eigen::VectorXd::Zero(robot.dimv())),
    beta(Eigen::VectorXd::Zero(robot.dimv())),
    mu_stack_(Eigen::VectorXd::Zero(robot.dim_passive()+robot.max_dimf())),
    f_stack_(Eigen::VectorXd::Zero(robot.max_dimf())),
    has_floating_base_(robot.has_floating_base()),
    dim_passive_(robot.dim_passive()),
    is_each_contact_active_(robot.max_point_contacts(), false),
    dimf_(robot.dimf()),
    dimc_(robot.dim_passive()+robot.dimf()) {
  robot.normalizeConfiguration(q);
}


inline SplitSolution::SplitSolution() 
  : lmd(),
    gmm(),
    mu_contact(),
    a(),
    f(),
    q(),
    v(),
    u(),
    beta(),
    mu_stack_(),
    f_stack_(),
    has_floating_base_(false),
    dim_passive_(0),
    is_each_contact_active_(),
    dimf_(0),
    dimc_(0) {
}


inline SplitSolution::~SplitSolution() {
}


inline void SplitSolution::setContactStatus(const Robot& robot) {
  is_each_contact_active_ = robot.is_each_contact_active();
  dimc_ = robot.dim_passive() + robot.dimf();
  dimf_ = robot.dimf();
}


inline Eigen::VectorBlock<Eigen::VectorXd> SplitSolution::mu_stack() {
  return mu_stack_.head(dimc_);
}


inline const Eigen::VectorBlock<const Eigen::VectorXd> 
SplitSolution::mu_stack() const {
  return mu_stack_.head(dimc_);
}


inline Eigen::VectorBlock<Eigen::VectorXd> SplitSolution::mu_floating_base() {
  return mu_stack_.head(dim_passive_);
}


inline const Eigen::VectorBlock<const Eigen::VectorXd> 
SplitSolution::mu_floating_base() const {
  return mu_stack_.head(dim_passive_);
}


inline Eigen::VectorBlock<Eigen::VectorXd> SplitSolution::mu_contacts() {
  return mu_stack_.segment(dim_passive_, dimf_);
}


inline const Eigen::VectorBlock<const Eigen::VectorXd> 
SplitSolution::mu_contacts() const {
  return mu_stack_.segment(dim_passive_, dimf_);
}


inline void SplitSolution::set_mu_stack() {
  int contact_index = 0;
  int segment_start = dim_passive_;
  for (const auto is_contact_active : is_each_contact_active_) {
    if (is_contact_active) {
      mu_stack_.template segment<3>(segment_start) = mu_contact[contact_index];
      segment_start += 3;
    }
    ++contact_index;
  }
}


inline void SplitSolution::set_mu_contacts() {
  int contact_index = 0;
  int segment_start = dim_passive_;
  for (const auto is_contact_active : is_each_contact_active_) {
    if (is_contact_active) {
      mu_contact[contact_index] = mu_stack_.template segment<3>(segment_start);
      segment_start += 3;
    }
    ++contact_index;
  }
}


inline void SplitSolution::set_f_stack() {
  int contact_index = 0;
  int segment_start = dim_passive_;
  for (const auto is_contact_active : is_each_contact_active_) {
    if (is_contact_active) {
      f_stack_.template segment<3>(segment_start) = f[contact_index];
      segment_start += 3;
    }
    ++contact_index;
  }
}


inline void SplitSolution::set_f() {
  int contact_index = 0;
  int segment_start = dim_passive_;
  for (const auto is_contact_active : is_each_contact_active_) {
    if (is_contact_active) {
      f[contact_index] = f_stack_.template segment<3>(segment_start);
      segment_start += 3;
    }
    ++contact_index;
  }
}


inline int SplitSolution::dimc() const {
  return dimc_;
}


inline int SplitSolution::dimf() const {
  return dimf_;
}

} // namespace idocp 

#endif // IDOCP_SPLIT_SOLUTION_HXX_