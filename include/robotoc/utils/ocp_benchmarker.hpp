#ifndef ROBOTOC_UTILS_OCP_BENCHMARKER_HPP_
#define ROBOTOC_UTILS_OCP_BENCHMARKER_HPP_ 

#include <memory>
#include <string>

#include "Eigen/Core"


namespace robotoc {
namespace benchmark {

template <typename OCPSolverType>
void CPUTime(OCPSolverType& ocp_solver, const double t, 
             const Eigen::VectorXd& q, const Eigen::VectorXd& v, 
             const int num_iteration=1000);

} // namespace benchmark
} // namespace robotoc 

#include "robotoc/utils/ocp_benchmarker.hxx"

#endif // ROBOTOC_UTILS_OCP_BENCHMARKER_HPP_