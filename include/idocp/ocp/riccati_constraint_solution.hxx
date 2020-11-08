#ifndef IDOCP_RICCATI_CONSTRAINT_SOLUTION_HXX_ 
#define IDOCP_RICCATI_CONSTRAINT_SOLUTION_HXX_

#include "idocp/ocp/riccati_constraint_solution.hpp"

#include <cassert>

namespace idocp {

inline RiccatiConstraintSolution::RiccatiConstraintSolution(
    const Robot& robot, const int N, const int max_impulse_stage) 
  : T_full_(N, Eigen::MatrixXd(2*robot.dimv(), robot.max_dimf())),
    T_impulse_full_(max_impulse_stage, 
                    Eigen::MatrixXd(2*robot.dimv(), robot.max_dimf())),
    EqNqq_full_(Eigen::MatrixXd(robot.max_dimf(), robot.dimv())),
    ENEt_full_(Eigen::MatrixXd(robot.max_dimf(), robot.max_dimf())),
    llt_(),
    dimv_(robot.dimv()),
    dimx_(2*robot.dimv()),
    dimf_(0),
    is_active_(false) { 
}


inline RiccatiConstraintSolution::RiccatiConstraintSolution() 
  : T_full_(),
    T_impulse_full_(),
    llt_(),
    EqNqq_full_(),
    ENEt_full_(),
    dimv_(0),
    dimx_(0),
    dimf_(0), 
    is_active_(false) { 
}


inline RiccatiConstraintSolution::~RiccatiConstraintSolution() { 
}


inline void RiccatiConstraintSolution::setImpulseStatus(
    const ImpulseStatus& impulse_status) {
  dimf_ = impulse_status.dimp();
  is_active_ = impulse_status.hasActiveImpulse();
}


inline bool RiccatiConstraintSolution::isActive() const {
  return is_active_;
}


inline Eigen::Block<Eigen::MatrixXd> RiccatiConstraintSolution::T(
    const int i) {
  return T_full_[i].topLeftCorner(dimx_, dimf_);
}


inline const Eigen::Block<const Eigen::MatrixXd> RiccatiConstraintSolution::T(
    const int i) const {
  return T_full_[i].topLeftCorner(dimx_, dimf_);
}


inline Eigen::Block<Eigen::MatrixXd> RiccatiConstraintSolution::T_impulse(
    const int i) {
  return T_impulse_full_[i].topLeftCorner(dimx_, dimf_);
}


inline const Eigen::Block<const Eigen::MatrixXd> 
RiccatiConstraintSolution::T_impulse(const int i) const {
  return T_impulse_full_[i].topLeftCorner(dimx_, dimf_);
}


inline Eigen::Block<Eigen::MatrixXd> RiccatiConstraintSolution::ENEt() {
  return ENEt_full_.topLeftCorner(dimf_, dimf_);
}


inline const Eigen::Block<const Eigen::MatrixXd> 
RiccatiConstraintSolution::ENEt() const {
  return ENEt_full_.topLeftCorner(dimf_, dimf_);
}


inline Eigen::Block<Eigen::MatrixXd> RiccatiConstraintSolution::EqNqq() {
  return EqNqq_full_.topLeftCorner(dimf_, dimv_);
}


inline const Eigen::Block<const Eigen::MatrixXd> 
RiccatiConstraintSolution::EqNqq() const {
  return EqNqq_full_.topLeftCorner(dimf_, dimv_);
}

} // namespace idocp

#endif // IDOCP_RICCATI_CONSTRAINT_SOLUTION_HXX_ 