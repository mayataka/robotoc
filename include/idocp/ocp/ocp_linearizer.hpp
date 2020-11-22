#ifndef IDOCP_OCP_LINEARIZER_HPP_ 
#define IDOCP_OCP_LINEARIZER_HPP_

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
  /// Must be non-negative. Default is 0.
  /// @param[in] num_proc Number of the threads in solving the optimal control 
  /// problem. Must be positive. Default is 1.
  ///
  OCPLinearizer(const double T, const int N, const int max_num_impulse=0, 
                const int num_proc=1);

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

  void initConstraints(HybridOCP& split_ocps, std::vector<Robot>& robots,
                       const ContactSequence& contact_sequence,
                       const double t, const HybridSolution& s) const;

  void linearizeOCP(HybridOCP& split_ocps, std::vector<Robot>& robots,
                    const ContactSequence& contact_sequence,
                    const double t, const Eigen::VectorXd& q, 
                    const Eigen::VectorXd& v, const HybridSolution& s,
                    HybridKKTMatrix& kkt_matrix,
                    HybridKKTResidual& kkt_residual) const;

  void computeKKTResidual(HybridOCP& split_ocps, std::vector<Robot>& robots, 
                          const ContactSequence& contact_sequence,
                          const double t, const Eigen::VectorXd& q, 
                          const Eigen::VectorXd& v, const HybridSolution& s,
                          HybridKKTMatrix& kkt_matrix, 
                          HybridKKTResidual& kkt_residual) const;

  double KKTError(const HybridOCP& split_ocps, 
                  const ContactSequence& contact_sequence, 
                  const HybridKKTResidual& kkt_residual);

  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

private:

  template <typename Algorithm>
  void runParallel(HybridOCP& split_ocps, std::vector<Robot>& robots,
                   const ContactSequence& contact_sequence,
                   const double t, const Eigen::VectorXd& q, 
                   const Eigen::VectorXd& v, const HybridSolution& s,
                   HybridKKTMatrix& kkt_matrix, 
                   HybridKKTResidual& kkt_residual) const;

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

#include "idocp/ocp/ocp_linearizer.hxx"

#endif // IDOCP_OCP_LINEARIZER_HPP_