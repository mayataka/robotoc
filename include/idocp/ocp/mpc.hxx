#ifndef IDOCP_MPC_HXX_
#define IDOCP_MPC_HXX_ 

#include "idocp/ocp/ocp.hpp"
// #include "idocp/ocp/parnmpc.hpp"

namespace idocp {

template <typename OCPType>
inline MPC<OCPType>::MPC(const Robot& robot, 
                         const std::shared_ptr<CostFunction>& cost,
                         const std::shared_ptr<Constraints>& constraints, 
                         const double T, const int N, const int num_proc)
  : ocp_(robot, cost, constraints, T, N, num_proc) {
}


template <typename OCPType>
inline MPC<OCPType>::MPC() 
  : ocp_() {
}


template <typename OCPType>
inline MPC<OCPType>::~MPC() {
}


template <>
inline void MPC<OCP>::initializeSolution(const double t, 
                                         const Eigen::VectorXd& q, 
                                         const Eigen::VectorXd& v, 
                                         const int max_itr) {
  ocp_.setStateTrajectory(q, v);
  for (int i=0; i<max_itr; ++i) {
    ocp_.updateSolution(t, q, v, true);
  }
}


// template <>
// inline void MPC<ParNMPC>::initializeSolution(const double t, 
//                                              const Eigen::VectorXd& q, 
//                                              const Eigen::VectorXd& v, 
//                                              const int max_itr) {
//   ocp_.setStateTrajectory(q, v);
//   ocp_.setAuxiliaryMatrixGuessByTerminalCost(t);
//   for (int i=0; i<max_itr; ++i) {
//     ocp_.updateSolution(t, q, v, true);
//   }
// }


template <typename OCPType>
inline void MPC<OCPType>::activateContact(const int contact_index, 
                                          const int time_stage_begin, 
                                          const int time_stage_end) {
  ocp_.activateContact(contact_index, time_stage_begin, time_stage_end);
}


template <typename OCPType>
inline void MPC<OCPType>::deactivateContact(const int contact_index, 
                                            const int time_stage_begin, 
                                            const int time_stage_end) {
  ocp_.deactivateContact(contact_index, time_stage_begin, time_stage_end);
}


template <typename OCPType>
inline void MPC<OCPType>::activateContacts(
    const std::vector<int>& contact_indices, const int time_stage_begin, 
    const int time_stage_end) {
  ocp_.activateContacts(contact_indices, time_stage_begin, time_stage_end);
}


template <typename OCPType>
inline void MPC<OCPType>::deactivateContacts(
    const std::vector<int>& contact_indices, const int time_stage_begin, 
    const int time_stage_end) {
  ocp_.deactivateContacts(contact_indices, time_stage_begin, time_stage_end);
}


template <typename OCPType>
inline void MPC<OCPType>::setContactPoint(
    const std::vector<Eigen::Vector3d>& contact_points) {
  ocp_.setContactSequence(contact_points);
}


template <typename OCPType>
inline void MPC<OCPType>::setContactPointByKinematics(
    const Eigen::VectorXd& q) {
  ocp_.setContactPointByKinematics(q);
}


template <typename OCPType>
inline void MPC<OCPType>::updateSolution(const double t, 
                                         const Eigen::VectorXd& q, 
                                         const Eigen::VectorXd& v) {
  ocp_.updateSolution(t, q, v, false);
}


template <typename OCPType>
inline void MPC<OCPType>::getControlInput(Eigen::VectorXd& u) {
  constexpr int initial_stage = 0;
  u.tail(12) = ocp_.getSolution(initial_stage).u;
}


template <typename OCPType>
inline void MPC<OCPType>::getStateFeedbackGain(Eigen::MatrixXd& Kq, 
                                               Eigen::MatrixXd& Kv) {
  constexpr int initial_stage = 0;
  ocp_.getStateFeedbackGain(initial_stage, Kq, Kv);
}


template <typename OCPType>
inline double MPC<OCPType>::KKTError() {
  return ocp_.KKTError();
}


template <typename OCPType>
inline void MPC<OCPType>::computeKKTResidual(const double t, 
                                             const Eigen::VectorXd& q, 
                                             const Eigen::VectorXd& v) {
  return ocp_.computeKKTResidual(t, q, v);
}

} // namespace idocp 

#endif // IDOCP_MPC_HXX_