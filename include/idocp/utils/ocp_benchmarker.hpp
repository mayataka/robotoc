#ifndef IDOCP_UTILS_OCP_BENCHMARKER_HPP_
#define IDOCP_UTILS_OCP_BENCHMARKER_HPP_ 

#include <memory>
#include <string>

#include "Eigen/Core"

#include "idocp/ocp/ocp.hpp"
#include "idocp/ocp/parnmpc.hpp"
#include "idocp/robot/robot.hpp"
#include "idocp/cost/cost_function.hpp"
#include "idocp/constraints/constraints.hpp"


namespace idocp {

template <typename OCPType>
class OCPBenchmarker {
public:
  OCPBenchmarker(const std::string& benchmark_name, const Robot& robot, 
                 const std::shared_ptr<CostFunction>& cost,
                 const std::shared_ptr<Constraints>& constraints, 
                 const double T, const int N, const int num_proc);

  OCPBenchmarker();

  ~OCPBenchmarker();

  // Use default copy constructor.
  OCPBenchmarker(const OCPBenchmarker&) = default;

  // Use default copy operator.
  OCPBenchmarker& operator=(const OCPBenchmarker&) = default;

  // Use default move constructor.
  OCPBenchmarker(OCPBenchmarker&&) noexcept = default;

  // Use default move operator.
  OCPBenchmarker& operator=(OCPBenchmarker&&) noexcept = default;

  void setInitialGuessSolution(const double t, const Eigen::VectorXd& q, 
                               const Eigen::VectorXd& v);

  ///
  /// @brief Activate a contact over specified time steps 
  /// (from start_time_stage to last_time_stage). 
  /// @param[in] contact_index Index of a contact of interedted. 
  /// @param[in] start_time_stage Start time stage. 
  /// @param[in] last_time_stage Last time stage. 
  ///
  void activateContact(const int contact_index, const int start_time_stage, 
                       const int last_time_stage);

  ///
  /// @brief Deactivate a contact over specified time steps 
  /// (from start_time_stage to last_time_stage). 
  /// @param[in] contact_index Index of a contact of interedted. 
  /// @param[in] start_time_stage Start time stage. 
  /// @param[in] last_time_stage Last time stage. 
  ///
  void deactivateContact(const int contact_index, const int start_time_stage, 
                         const int last_time_stage);

  ///
  /// @brief Activate contacts over specified time steps 
  /// (from start_time_stage to last_time_stage). 
  /// @param[in] contact_indices Indices of contacts of interedted. 
  /// @param[in] start_time_stage Start time stage. 
  /// @param[in] last_time_stage Last time stage. 
  ///
  void activateContacts(const std::vector<int>& contact_indices, 
                        const int start_time_stage, const int last_time_stage);

  ///
  /// @brief Deactivate contacts over specified time steps 
  /// (from start_time_stage to last_time_stage). 
  /// @param[in] contact_indices Indices of contacts of interedted. 
  /// @param[in] start_time_stage Start time stage. 
  /// @param[in] last_time_stage Last time stage. 
  ///
  void deactivateContacts(const std::vector<int>& contact_indices, 
                          const int start_time_stage, const int last_time_stage);

  void testCPUTime(const double t, const Eigen::VectorXd& q, 
                   const Eigen::VectorXd& v, const int num_iteration=1000,
                   const bool line_search=false);

  void testConvergence(const double t, const Eigen::VectorXd& q, 
                       const Eigen::VectorXd& v, const int num_iteration=10,
                       const bool line_search=false);

  void printSolution();

private:
  std::string benchmark_name_;
  OCPType ocp_;
  int dimq_, dimv_, max_dimf_, N_, num_proc_;
  double T_, dtau_;

};

} // namespace idocp 

#include "ocp_benchmarker.hxx"

#endif // IDOCP_UTILS_OCP_BENCHMARKER_HPP_