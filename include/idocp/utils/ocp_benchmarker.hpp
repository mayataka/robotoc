#ifndef IDOCP_UTILS_OCP_BENCHMARKER_HPP_
#define IDOCP_UTILS_OCP_BENCHMARKER_HPP_ 

#include <memory>
#include <string>

#include "Eigen/Core"

#include "idocp/utils/logger.hpp"


namespace idocp {
namespace benchmark {

template <typename OCPSolverType>
void CPUTime(OCPSolverType& ocp_solver, const double t, 
             const Eigen::VectorXd& q, const Eigen::VectorXd& v, 
             const int num_iteration=1000, const bool line_search=false);

template <typename OCPSolverType>
void Convergence(OCPSolverType& ocp_solver, const double t, 
                 const Eigen::VectorXd& q, const Eigen::VectorXd& v, 
                 const int num_iteration=10, const bool line_search=false);

template <typename OCPSolverType>
void Convergence(OCPSolverType& ocp_solver, Logger& logger, 
                 const double t, const Eigen::VectorXd& q, 
                 const Eigen::VectorXd& v, const int num_iteration=10, 
                 const bool line_search=false);

} // namespace benchmark
} // namespace idocp 

#include "idocp/utils/ocp_benchmarker.hxx"

#endif // IDOCP_UTILS_OCP_BENCHMARKER_HPP_