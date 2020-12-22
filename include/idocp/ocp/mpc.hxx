#ifndef IDOCP_MPC_HXX_
#define IDOCP_MPC_HXX_ 

#include "idocp/ocp/ocp_solver.hpp"
// #include "idocp/ocp/parnmpc_solver.hpp"

namespace idocp {

template <typename OCPSolverType>
inline MPC<OCPSolverType>::MPC(
    const Robot& robot, const std::shared_ptr<CostFunction>& cost, 
    const std::shared_ptr<Constraints>& constraints, const double T, 
    const int N, const int max_num_impulse, const int num_proc)
  : ocp_solver_(robot, cost, constraints, T, N, max_num_impulse, num_proc) {
}


template <typename OCPSolverType>
inline MPC<OCPSolverType>::MPC() 
  : ocp_solver_() {
}


template <typename OCPSolverType>
inline MPC<OCPSolverType>::~MPC() {
}


template <>
inline void MPC<OCPSolver>::initializeSolution(const double t, 
                                               const Eigen::VectorXd& q, 
                                               const Eigen::VectorXd& v, 
                                               const int max_itr) {
  ocp_solver_.setStateTrajectory(t, q, v);
  for (int i=0; i<max_itr; ++i) {
    ocp_solver_.updateSolution(t, q, v, true);
  }
}


// template <>
// inline void MPC<ParNMPCSolver>::initializeSolution(const double t, 
//                                                    const Eigen::VectorXd& q, 
//                                                    const Eigen::VectorXd& v, 
//                                                    const int max_itr) {
//   ocp_solver_.setStateTrajectory(t, q, v);
//   ocp_solver_.setAuxiliaryMatrixGuessByTerminalCost(t);
//   for (int i=0; i<max_itr; ++i) {
//     ocp_solver_.updateSolution(t, q, v, true);
//   }
// }


template <typename OCPSolverType>
inline void MPC<OCPSolverType>::updateSolution(
    const double t, const Eigen::VectorXd& q, const Eigen::VectorXd& v,
    const int max_iter, const double KKT_tol) {
  for (int i=0; i<max_iter; ++i) {
    ocp_solver_.updateSolution(t, q, v, false);
    if (KKT_tol > 0) {
      ocp_solver_.computeKKTResidual(t, q, v);
      if (ocp_solver_.KKTError() < KKT_tol) {
        break;
      }
    }
  }
}


template <typename OCPSolverType>
inline void MPC<OCPSolverType>::updateSolutionWithContinuationMethod(
    const double t, const Eigen::VectorXd& q, const Eigen::VectorXd& v, 
    const double sampling_period) {
  assert(sampling_period > 0);
  ocp_solver_.updateSolutionWithContinuationMethod(t, q, v, sampling_period);
}


template <typename OCPSolverType>
inline void MPC<OCPSolverType>::getControlInput(Eigen::VectorXd& u) {
  u = ocp_solver_.getSolution(0).u;
}


template <typename OCPSolverType>
inline void MPC<OCPSolverType>::getStateFeedbackGain(Eigen::MatrixXd& Kq, 
                                                     Eigen::MatrixXd& Kv) {
  constexpr int initial_stage = 0;
  ocp_solver_.getStateFeedbackGain(initial_stage, Kq, Kv);
}


template <typename OCPSolverType>
inline double MPC<OCPSolverType>::KKTError() {
  return ocp_solver_.KKTError();
}


template <typename OCPSolverType>
inline void MPC<OCPSolverType>::computeKKTResidual(const double t, 
                                             const Eigen::VectorXd& q, 
                                             const Eigen::VectorXd& v) {
  return ocp_solver_.computeKKTResidual(t, q, v);
}


template <typename OCPSolverType>
inline OCPSolverType* MPC<OCPSolverType>::getSolverHandle() {
  return &ocp_solver_;
}

} // namespace idocp 

#endif // IDOCP_MPC_HXX_