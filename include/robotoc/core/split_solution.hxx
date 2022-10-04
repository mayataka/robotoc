#ifndef ROBOTOC_SPLIT_SOLUTION_HXX_
#define ROBOTOC_SPLIT_SOLUTION_HXX_

#include "robotoc/core/split_solution.hpp"

#include <cassert>

namespace robotoc {

inline void SplitSolution::setContactStatus(
    const ContactStatus& contact_status) {
  assert(contact_status.maxNumContacts() == is_contact_active_.size());
  is_contact_active_ = contact_status.isContactActive();
  dimf_ = contact_status.dimf();
}


inline void SplitSolution::setContactStatus(
    const ImpactStatus& contact_status) {
  assert(contact_status.maxNumContacts() == is_contact_active_.size());
  is_contact_active_ = contact_status.isImpactActive();
  dimf_ = contact_status.dimf();
}


inline void SplitSolution::setContactStatus(const SplitSolution& other) {
  assert(other.isContactActive().size() == is_contact_active_.size());
  is_contact_active_ = other.isContactActive();
  dimf_ = other.dimf();
}


inline void SplitSolution::setSwitchingConstraintDimension(const int dims) {
  assert(dims >= 0);
  assert(dims <= xi_stack_.size());
  dims_ = std::max(dims, 0);
}


inline Eigen::VectorBlock<Eigen::VectorXd> SplitSolution::f_stack() {
  return f_stack_.head(dimf_);
}


inline const Eigen::VectorBlock<const Eigen::VectorXd> 
SplitSolution::f_stack() const {
  return f_stack_.head(dimf_);
}


inline void SplitSolution::set_f_stack() {
  int segment_begin = 0;
  for (int i=0; i<max_num_contacts_; ++i) {
    if (is_contact_active_[i]) {
      switch (contact_types_[i]) {
      case ContactType::PointContact:
        f_stack_.template segment<3>(segment_begin) = f[i].template head<3>();
        segment_begin += 3;
        break;
      case ContactType::SurfaceContact:
        f_stack_.template segment<6>(segment_begin) = f[i];
        segment_begin += 6;
        break;
      default:
        break;
      }
    }
  }
}


inline void SplitSolution::set_f_vector() {
  int segment_begin = 0;
  for (int i=0; i<max_num_contacts_; ++i) {
    if (is_contact_active_[i]) {
      switch (contact_types_[i]) {
      case ContactType::PointContact:
        f[i].template head<3>() = f_stack_.template segment<3>(segment_begin);
        segment_begin += 3;
        break;
      case ContactType::SurfaceContact:
        f[i] = f_stack_.template segment<6>(segment_begin);
        segment_begin += 6;
        break;
      default:
        break;
      }
    }
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
  int segment_begin = 0;
  for (int i=0; i<max_num_contacts_; ++i) {
    if (is_contact_active_[i]) {
      switch (contact_types_[i]) {
      case ContactType::PointContact:
        mu_stack_.template segment<3>(segment_begin) = mu[i].template head<3>();
        segment_begin += 3;
        break;
      case ContactType::SurfaceContact:
        mu_stack_.template segment<6>(segment_begin) = mu[i];
        segment_begin += 6;
        break;
      default:
        break;
      }
    }
  }
}


inline void SplitSolution::set_mu_vector() {
  int segment_begin = 0;
  for (int i=0; i<max_num_contacts_; ++i) {
    if (is_contact_active_[i]) {
      switch (contact_types_[i]) {
      case ContactType::PointContact:
        mu[i].template head<3>() = mu_stack_.template segment<3>(segment_begin);
        segment_begin += 3;
        break;
      case ContactType::SurfaceContact:
        mu[i] = mu_stack_.template segment<6>(segment_begin);
        segment_begin += 6;
        break;
      default:
        break;
      }
    }
  }
}


inline Eigen::VectorBlock<Eigen::VectorXd> SplitSolution::xi_stack() {
  return xi_stack_.head(dims_);
}


inline const Eigen::VectorBlock<const Eigen::VectorXd> 
SplitSolution::xi_stack() const {
  return xi_stack_.head(dims_);
}


inline int SplitSolution::dimf() const {
  return dimf_;
}


inline int SplitSolution::dims() const {
  return dims_;
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

} // namespace robotoc 

#endif // ROBOTOC_SPLIT_SOLUTION_HXX_