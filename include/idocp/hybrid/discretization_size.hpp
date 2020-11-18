#ifndef IDOCP_DISCRETIZATION_SIZE_HPP_
#define IDOCP_DISCRETIZATION_SIZE_HPP_

#include <vector>
#include <cassert>

#include "idocp/hybrid/contact_sequence.hpp"

namespace idocp {

class DiscretizationSize {
public:
  DiscretizationSize(const double T, const int N) 
    : T_(T),
      dtau_origin_(T/N),
      N_(N),
      total_num_impulse_(0),
      total_num_lift_(0),
      dtau_(N, T/N),
      dtau_aux_(N, 0),
      dtau_lift_(N, 0) {
  }

  DiscretizationSize() 
    : T_(0),
      dtau_origin_(0),
      N_(0),
      total_num_impulse_(0),
      total_num_lift_(0),
      dtau_(),
      dtau_aux_(),
      dtau_lift_() {
  }

  ~DiscretizationSize() {
  }

  void setContactSequence(const ContactSequence& contact_sequence) {
    const double dtau = T_ / N_;
    for (int i=0; i<N_; ++i) {
      if (contact_sequence.existImpulseStage(i)) {
        dtau_[i] = contact_sequence.impulseTime(contact_sequence.impulseIndex(i)) 
                    - i * dtau_origin_;
        dtau_aux_[i] = dtau_origin_ - dtau_[i];
      }
      else if (contact_sequence.existLiftStage(i)) {
        dtau_[i] = contact_sequence.liftTime(contact_sequence.liftIndex(i)) 
                    - i * dtau_origin_;
        dtau_lift_[i] = dtau_origin_ - dtau_[i];
      }
      else {
        dtau_[i] = dtau;
      }
    }
  }

  double operator[](const int time_stage) const {
    assert(time_stage >= 0);
    assert(time_stage < N_);
    return dtau_[time_stage];
  }

  double aux(const int impulse_index) const {
    assert(impulse_index >= 0);
    assert(impulse_index < N_);
    return dtau_aux_[impulse_index];
  }

  double lift(const int lift_index) const {
    assert(lift_index >= 0);
    assert(lift_index < N_);
    return dtau_lift_[lift_index];
  }

private:
  double T_, dtau_origin_;
  int N_, total_num_impulse_, total_num_lift_;
  std::vector<double> dtau_;
  std::vector<double> dtau_aux_;
  std::vector<double> dtau_lift_;

};
  
} // namespace idocp

#endif // IDOCP_DISCRETIZATION_SIZE_HPP_ 