#ifndef IDOCP_MPC_HXX_
#define IDOCP_MPC_HXX_ 

#include "idocp/ocp/ocp.hpp"
// #include "idocp/ocp/parnmpc.hpp"

namespace idocp {

template <typename OCPType>
inline MPC<OCPType>::MPC(const Robot& robot, 
                         const std::shared_ptr<CostFunction>& cost,
                         const std::shared_ptr<Constraints>& constraints, 
                         const double T, const int N, const int max_num_impulse,
                         const int num_proc)
  : ocp_(robot, cost, constraints, T, N, max_num_impulse, num_proc) {
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
inline void MPC<OCPType>::updateSolution(const double t, 
                                         const Eigen::VectorXd& q, 
                                         const Eigen::VectorXd& v,
                                         const int max_iter,
                                         const double KKT_tol) {
  for (int i=0; i<max_iter; ++i) {
    ocp_.updateSolution(t, q, v, false);
    if (KKT_tol > 0) {
      ocp_.computeKKTResidual(t, q, v);
      if (ocp_.KKTError() < KKT_tol) {
        break;
      }
    }
  }
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


template <typename OCPType>
inline OCPType* MPC<OCPType>::getSolverHandle() {
  return &ocp_;
}

} // namespace idocp 

#endif // IDOCP_MPC_HXX_