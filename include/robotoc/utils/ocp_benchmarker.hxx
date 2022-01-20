#ifndef ROBOTOC_UTILS_OCP_BENCHMARKER_HXX_
#define ROBOTOC_UTILS_OCP_BENCHMARKER_HXX_ 

#include "robotoc/utils/ocp_benchmarker.hxx"

#include <iostream>
#include <chrono>


namespace robotoc {
namespace benchmark {

template <typename OCPSolverType>
inline void CPUTime(OCPSolverType& ocp_solver, const double t, 
                    const Eigen::VectorXd& q, const Eigen::VectorXd& v, 
                    const int num_iteration) {
  std::chrono::high_resolution_clock::time_point start_clock, end_clock;
  start_clock = std::chrono::high_resolution_clock::now();
  for (int i=0; i<num_iteration; ++i) {
    ocp_solver.updateSolution(t, q, v);
  }
  end_clock = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double, std::milli> timing = end_clock - start_clock;
  std::cout << "---------- OCP benchmark : CPU time ----------" << std::endl;
  std::cout << "total CPU time: " << timing.count() 
            << "[ms]" << std::endl;
  std::cout << "CPU time per update: " << timing.count() / num_iteration
            << "[ms]" << std::endl;
  std::cout << "-----------------------------------" << std::endl;
  std::cout << std::endl;
}

} // namespace benchmark
} // namespace robotoc 

#endif // ROBOTOC_UTILS_OCP_BENCHMARKER_HXX_