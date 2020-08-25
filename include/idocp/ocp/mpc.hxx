#ifndef IDOCP_MPC_HXX_
#define IDOCP_MPC_HXX_ 

#include "idocp/ocp/ocp.hpp"
#include "idocp/ocp/parnmpc.hpp"

namespace idocp {

template <typename OCPType>
inline MPC<OCPType>::MPC(const Robot& robot, 
                         const std::shared_ptr<CostFunction>& cost,
                         const std::shared_ptr<Constraints>& constraints, 
                         const double T, const int N, const int num_proc)
  : ocp_(robot, cost, constraints, T, N, num_proc) {
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


template <>
inline void MPC<ParNMPC>::initializeSolution(const double t, 
                                             const Eigen::VectorXd& q, 
                                             const Eigen::VectorXd& v, 
                                             const int max_itr) {
  ocp_.setStateTrajectory(q, v);
  ocp_.setAuxiliaryMatrixGuessByTerminalCost(t);
  for (int i=0; i<max_itr; ++i) {
    ocp_.updateSolution(t, q, v, true);
  }
}


template <typename OCPType>
inline void MPC<OCPType>::updateSolution(const double t, 
                                         const Eigen::VectorXd& q, 
                                         const Eigen::VectorXd& v) {
  ocp_.updateSolution(t, q, v, false);
}


template <typename OCPType>
inline void MPC<OCPType>::getControlInput(Eigen::VectorXd& u) {
  ocp_.getInitialControlInput(u);
}


template <typename OCPType>
inline void MPC<OCPType>::getStateFeedbackGain(Eigen::MatrixXd& Kq, 
                                               Eigen::MatrixXd& Kv) {
  ocp_.getStateFeedbackGain(Kq, Kv);
}


template <typename OCPType>
inline double MPC<OCPType>::KKTError(const double t, const Eigen::VectorXd& q, 
                                     const Eigen::VectorXd& v) {
  return ocp_.KKTError(t, q, v);
}

} // namespace idocp 

#endif // IDOCP_MPC_HXX_