#ifndef IDOCP_OCP_LINEARIZER_HPP_ 
#define IDOCP_OCP_LINEARIZER_HPP_

#include <vector>

#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/hybrid/hybrid_container.hpp"
#include "idocp/hybrid/contact_sequence.hpp"


namespace idocp {

///
/// @class OCPLinearizer
/// @brief Linearize of the optimal control problem. 
///
class OCPLinearizer {
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
  OCPLinearizer(const double T, const int N, const int max_num_impulse, 
                const int num_proc);

  ///
  /// @brief Default constructor. 
  ///
  OCPLinearizer();

  ///
  /// @brief Destructor. 
  ///
  ~OCPLinearizer();

  ///
  /// @brief Default copy constructor. 
  ///
  OCPLinearizer(const OCPLinearizer&) = default;

  ///
  /// @brief Default copy assign operator. 
  ///
  OCPLinearizer& operator=(const OCPLinearizer&) = default;

  ///
  /// @brief Default move constructor. 
  ///
  OCPLinearizer(OCPLinearizer&&) noexcept = default;

  ///
  /// @brief Default move assign operator. 
  ///
  OCPLinearizer& operator=(OCPLinearizer&&) noexcept = default;

  void initConstraints(OCP& ocp, std::vector<Robot>& robots,
                       const ContactSequence& contact_sequence, 
                       const Solution& s) const;

  void linearizeOCP(OCP& ocp, std::vector<Robot>& robots,
                    const ContactSequence& contact_sequence,
                    const double t, const Eigen::VectorXd& q, 
                    const Eigen::VectorXd& v, const Solution& s,
                    KKTMatrix& kkt_matrix, KKTResidual& kkt_residual) const;

  void computeKKTResidual(OCP& ocp, std::vector<Robot>& robots, 
                          const ContactSequence& contact_sequence,
                          const double t, const Eigen::VectorXd& q, 
                          const Eigen::VectorXd& v, const Solution& s,
                          KKTMatrix& kkt_matrix, 
                          KKTResidual& kkt_residual) const;

  double KKTError(const OCP& ocp, const ContactSequence& contact_sequence, 
                  const KKTResidual& kkt_residual);

  void integrateSolution(OCP& ocp, const std::vector<Robot>& robots,
                         const ContactSequence& contact_sequence,
                         const KKTMatrix& kkt_matrix,
                         const KKTResidual& kkt_residual,
                         const double primal_step_size,
                         const double dual_step_size,
                         Direction& d, Solution& s) const;

  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

private:

  template <typename Algorithm>
  void runParallel(OCP& ocp, std::vector<Robot>& robots,
                   const ContactSequence& contact_sequence,
                   const double t, const Eigen::VectorXd& q, 
                   const Eigen::VectorXd& v, const Solution& s,
                   KKTMatrix& kkt_matrix, KKTResidual& kkt_residual) const;

  const Eigen::VectorXd& q_prev(const ContactSequence& contact_sequence, 
                                const Eigen::VectorXd& q,
                                const Solution& s, const int time_stage) const;

  double dtau(const ContactSequence& contact_sequence, 
              const int time_stage) const;

  static bool is_state_constraint_valid(const int time_stage_before_impulse);

  double T_, dtau_;
  int N_, num_proc_;
  Eigen::VectorXd kkt_error_;
};

} // namespace idocp 

#include "idocp/ocp/ocp_linearizer.hxx"

#endif // IDOCP_OCP_LINEARIZER_HPP_