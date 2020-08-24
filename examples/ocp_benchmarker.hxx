#ifndef IDOCP_OCP_BENCHMARKER_HXX_
#define IDOCP_OCP_BENCHMARKER_HXX_ 

#include <iostream>
#include <chrono>

#include "idocp/ocp/ocp.hpp"
#include "idocp/ocp/parnmpc.hpp"

namespace idocp {

template <typename OCPType>
inline OCPBenchmarker<OCPType>::OCPBenchmarker(
    const std::string& benchmark_name, const Robot& robot, 
    const std::shared_ptr<CostFunction>& cost, 
    const std::shared_ptr<Constraints>& constraints, const double T, 
    const int N, const int num_proc)
  : benchmark_name_(benchmark_name),
    ocp_(robot, cost, constraints, T, N, num_proc),
    dimq_(robot.dimq()),
    dimv_(robot.dimv()),
    max_dimf_(robot.max_dimf()),
    N_(N),
    num_proc_(num_proc),
    T_(T),
    dtau_(T/N) {
}


template <typename OCPType>
inline OCPBenchmarker<OCPType>::OCPBenchmarker()
  : benchmark_name_(),
    ocp_(),
    dimq_(0),
    dimv_(0),
    max_dimf_(0),
    N_(0),
    num_proc_(0),
    T_(0),
    dtau_(0) {
}


template <typename OCPType>
inline OCPBenchmarker<OCPType>::~OCPBenchmarker() {
}


template <>
inline void OCPBenchmarker<OCP>::setInitialGuessSolution(
    const double t, const Eigen::VectorXd& q, const Eigen::VectorXd& v) {
  ocp_.setStateTrajectory(q, v);
}


template <>
inline void OCPBenchmarker<ParNMPC>::setInitialGuessSolution(
    const double t, const Eigen::VectorXd& q, const Eigen::VectorXd& v) {
  ocp_.setStateTrajectory(q, v);
  ocp_.setAuxiliaryMatrixGuessByTerminalCost(t);
}


template <typename OCPType>
inline void OCPBenchmarker<OCPType>::setContactStatus(
    const std::vector<bool>& contact_status) {
  const auto contact_sequence = std::vector<std::vector<bool>>(N_, contact_status);
  ocp_.setContactSequence(contact_sequence);
}


template <typename OCPType>
inline void OCPBenchmarker<OCPType>::testCPUTime(const double t, 
                                                 const Eigen::VectorXd& q, 
                                                 const Eigen::VectorXd& v,
                                                 const int num_iteration,
                                                 const bool line_search) {
  std::chrono::system_clock::time_point start_clock, end_clock;
  start_clock = std::chrono::system_clock::now();
  for (int i=0; i<num_iteration; ++i) {
    ocp_.updateSolution(t, q, v, line_search);
  }
  end_clock = std::chrono::system_clock::now();
  std::cout << "---------- OCP benchmark ----------" << std::endl;
  std::cout << "Test CPU time of " << benchmark_name_ << std::endl;
  std::cout << "robot.dimq() = " << dimq_ << std::endl;
  std::cout << "robot.dimv() = " << dimv_ << std::endl;
  std::cout << "robot.max_dimf() = " << max_dimf_ << std::endl;
  std::cout << "N (number of stages) = " << N_ << std::endl;
  if (line_search) {
    std::cout << "Line search is enable" << std::endl;
  }
  else {
    std::cout << "Line search is disable" << std::endl;
  }
  std::cout << "number of threads = " << num_proc_ << std::endl;
  std::cout << "total CPU time: " 
            << 1e-03 * std::chrono::duration_cast<std::chrono::microseconds>(
                  end_clock-start_clock).count() 
            << "[ms]" << std::endl;
  std::cout << "CPU time per update: " 
            << 1e-03 * std::chrono::duration_cast<std::chrono::microseconds>(
                  end_clock-start_clock).count() / num_iteration 
            << "[ms]" << std::endl;
  std::cout << "-----------------------------------" << std::endl;
  std::cout << std::endl;
}


template <typename OCPType>
inline void OCPBenchmarker<OCPType>::testConvergence(const double t, 
                                                     const Eigen::VectorXd& q, 
                                                     const Eigen::VectorXd& v,
                                                     const int num_iteration,
                                                     const bool line_search) {
  std::cout << "---------- OCP benchmark ----------" << std::endl;
  std::cout << "Test convergence of " << benchmark_name_ << std::endl;
  std::cout << "robot.dimq() = " << dimq_ << std::endl;
  std::cout << "robot.dimv() = " << dimv_ << std::endl;
  std::cout << "robot.max_dimf() = " << max_dimf_ << std::endl;
  std::cout << "N (number of stages) = " << N_ << std::endl;
  std::cout << "q = " << q.transpose() << std::endl;
  std::cout << "v = " << v.transpose() << std::endl;
  if (line_search) {
    std::cout << "Line search is enable" << std::endl;
  }
  else {
    std::cout << "Line search is disable" << std::endl;
  }
  std::cout << "Initial KKT error = " << ocp_.KKTError(t, q, v) << std::endl;
  for (int i=0; i<num_iteration; ++i) {
    ocp_.updateSolution(t, q, v, line_search);
    std::cout << "KKT error at iteration " << i << " = " 
              << ocp_.KKTError(t, q, v) << std::endl;
  }
  std::cout << "-----------------------------------" << std::endl;
  std::cout << std::endl;
}


template <typename OCPType>
inline void OCPBenchmarker<OCPType>::printSolution() {
  ocp_.printSolution();
}

} // namespace idocp 

#endif // IDOCP_OCP_BENCHMARKER_HXX_