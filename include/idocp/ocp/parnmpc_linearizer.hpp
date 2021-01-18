#ifndef IDOCP_PARNMPC_LINEARIZER_HPP_
#define IDOCP_PARNMPC_LINEARIZER_HPP_

#include <vector>
#include <memory>

#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/hybrid/hybrid_container.hpp"
#include "idocp/hybrid/contact_sequence.hpp"


namespace idocp {

///
/// @class ParNMPCLinearizer
/// @brief Linearize of the optimal control problem. 
///
class ParNMPCLinearizer {
public:
  ///
  /// @brief Construct optimal control problem solver.
  /// @param[in] N Number of discretization of the horizon. Must be more than 1. 
  /// @param[in] max_num_impulse Maximum number of the impulse on the horizon. 
  /// Must be non-negative. 
  /// @param[in] nthreads Number of the threads in solving the optimal control 
  /// problem. Must be positive. 
  ///
  ParNMPCLinearizer(const int N, const int max_num_impulse, const int nthreads);

  ///
  /// @brief Default constructor. 
  ///
  ParNMPCLinearizer();

  ///
  /// @brief Destructor. 
  ///
  ~ParNMPCLinearizer();

  ///
  /// @brief Default copy constructor. 
  ///
  ParNMPCLinearizer(const ParNMPCLinearizer&) = default;

  ///
  /// @brief Default copy assign operator. 
  ///
  ParNMPCLinearizer& operator=(const ParNMPCLinearizer&) = default;

  ///
  /// @brief Default move constructor. 
  ///
  ParNMPCLinearizer(ParNMPCLinearizer&&) noexcept = default;

  ///
  /// @brief Default move assign operator. 
  ///
  ParNMPCLinearizer& operator=(ParNMPCLinearizer&&) noexcept = default;

  void initConstraints(ParNMPC& parnmpc, std::vector<Robot>& robots,
                       const ContactSequence& contact_sequence, 
                       const Solution& s) const;

  void computeKKTResidual(ParNMPC& parnmpc, std::vector<Robot>& robots, 
                          const ContactSequence& contact_sequence,
                          const Eigen::VectorXd& q, const Eigen::VectorXd& v, 
                          const Solution& s, KKTMatrix& kkt_matrix, 
                          KKTResidual& kkt_residual) const;

  double KKTError(const ParNMPC& parnmpc, const KKTResidual& kkt_residual);

  void printKKTError() const;

  void integrateSolution(ParNMPC& parnmpc, const std::vector<Robot>& robots,
                         const KKTMatrix& kkt_matrix,
                         const KKTResidual& kkt_residual,
                         const double primal_step_size,
                         const double dual_step_size,
                         const Direction& d, Solution& s) const;

  static const Eigen::VectorXd& q_prev(const OCPDiscretizer& ocp_discretizer, 
                                       const Eigen::VectorXd& q, 
                                       const Solution& s, const int time_stage);

  static const Eigen::VectorXd& v_prev(const OCPDiscretizer& ocp_discretizer, 
                                       const Eigen::VectorXd& v, 
                                       const Solution& s, const int time_stage);


  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

private:

  static constexpr double kMindtau
      = std::sqrt(std::numeric_limits<double>::epsilon());

  int N_, max_num_impulse_, nthreads_;
  Eigen::VectorXd kkt_error_;

};

} // namespace idocp 

#include "idocp/ocp/parnmpc_linearizer.hxx"

#endif // IDOCP_PARNMPC_LINEARIZER_HPP_ 