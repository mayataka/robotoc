#ifndef IDOCP_OCP_PARALLEL_BASE_HPP_ 
#define IDOCP_OCP_PARALLEL_BASE_HPP_

#include <vector>

#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/ocp/split_ocp.hpp"
#include "idocp/impulse/split_impulse_ocp.hpp"
#include "idocp/ocp/terminal_ocp.hpp"
#include "idocp/ocp/split_solution.hpp"
#include "idocp/ocp/kkt_matrix.hpp"
#include "idocp/ocp/kkt_residual.hpp"
#include "idocp/impulse/impulse_split_solution.hpp"
#include "idocp/impulse/impulse_kkt_matrix.hpp"
#include "idocp/impulse/impulse_kkt_residual.hpp"
#include "idocp/hybrid/hybrid_container.hpp"
#include "idocp/hybrid/contact_sequence.hpp"


namespace idocp {

///
/// @class OCPParallelBase
/// @brief Linearize of the optimal control problem. 
///
template <typename Derived>
class OCPParallelBase {
public:
  ///
  /// @brief Construct optimal control problem solver.
  /// @param[in] T Length of the horizon. Must be positive.
  /// @param[in] N Number of discretization of the horizon. Must be more than 1. 
  /// @param[in] max_num_impulse Maximum number of the impulse on the horizon. 
  /// Must be non-negative. 
  /// @param[in] num_proc Number of the threads in solving the optimal control 
  /// problem. Must be positive. 
  ///
  OCPParallelBase(const double T, const int N, const int max_num_impulse, 
                  const int num_proc);

  ///
  /// @brief Default constructor. 
  ///
  OCPParallelBase();

  ///
  /// @brief Destructor. 
  ///
  ~OCPParallelBase();

  ///
  /// @brief Default copy constructor. 
  ///
  OCPParallelBase(const OCPParallelBase&) = default;

  ///
  /// @brief Default copy assign operator. 
  ///
  OCPParallelBase& operator=(const OCPParallelBase&) = default;

  ///
  /// @brief Default move constructor. 
  ///
  OCPParallelBase(OCPParallelBase&&) noexcept = default;

  ///
  /// @brief Default move assign operator. 
  ///
  OCPParallelBase& operator=(OCPParallelBase&&) noexcept = default;

  template <typename... Args>
  void runParallel(Args... args);

  template <typename... Args>
  void runParallel(Args... args) const;

  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

private:

  const Eigen::VectorXd& q_prev(const ContactSequence& contact_sequence, 
                                const Eigen::VectorXd& q,
                                const HybridSolution& s,
                                const int time_stage) const;

  double dtau(const ContactSequence& contact_sequence, 
              const int time_stage) const;


  double T_, dtau_;
  int N_, num_proc_;
  Eigen::VectorXd kkt_error_;
};

} // namespace idocp 

#include "idocp/ocp/ocp_parallel_base.hxx"

#endif // IDOCP_OCP_PARALLEL_BASE_HPP_