#ifndef IDOCP_OCP_LINEARIZER_HPP_ 
#define IDOCP_OCP_LINEARIZER_HPP_

#include <vector>
#include <memory>

#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/cost/cost_function.hpp"
#include "idocp/constraints/constraints.hpp"
#include "idocp/cost/impulse_cost_function.hpp"
#include "idocp/constraints/impulse_constraints.hpp"
#include "idocp/ocp/split_ocp.hpp"
#include "idocp/impulse/split_impulse_ocp.hpp"
#include "idocp/ocp/terminal_ocp.hpp"
#include "idocp/ocp/split_solution.hpp"
#include "idocp/ocp/split_direction.hpp"
#include "idocp/impulse/impulse_split_solution.hpp"
#include "idocp/impulse/impulse_split_direction.hpp"
#include "idocp/hybrid/hybrid_container.hpp"
#include "idocp/hybrid/contact_sequence.hpp"


namespace idocp {

///
/// @class OCPLinearizer
/// @brief Linearize of the optimal control problem. 
///
class OCPLinearizer {
public:
  using HybridOCP = hybrid_container<SplitOCP, SplitImpulseOCP>;
  using HybridSolution = hybrid_container<SplitSolution, ImpulseSplitSolution>;
  using HybridDirection = hybrid_container<SplitDirection, ImpulseSplitDirection>;
  using HybridKKTMatrix = hybrid_container<KKTMatrix, ImpulseKKTMatrix>;
  using HybridKKTResidual = hybrid_container<KKTResidual, ImpulseKKTResidual>;

  ///
  /// @brief Construct optimal control problem solver.
  /// @param[in] robot Robot model. Must be initialized by URDF or XML.
  /// @param[in] cost Shared ptr to the cost function.
  /// @param[in] constraints Shared ptr to the constraints.
  /// @param[in] T Length of the horizon. Must be positive.
  /// @param[in] N Number of discretization of the horizon. Must be more than 1. 
  /// @param[in] max_num_impulse Maximum number of the impulse on the horizon. 
  /// Must be non-negative. Default is 0.
  /// @param[in] num_proc Number of the threads in solving the optimal control 
  /// problem. Must be positive. Default is 1.
  ///
  OCPLinearizer(const Robot& robot, const std::shared_ptr<CostFunction>& cost,
                const std::shared_ptr<Constraints>& constraints, const double T, 
                const int N, const int max_num_impulse=0, const int num_proc=1);

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

  void linearizeOCP(std::vector<Robot>& robots,
                    const ContactSequence& contact_sequence,
                    const double t, const Eigen::VectorXd& q, 
                    const Eigen::VectorXd& v, const HybridSolution& s,
                    HybridKKTMatrix& kkt_matrix,
                    HybridKKTResidual& kkt_residual);

  void linearizeInitialState(const Robot& robot, const double t, 
                             const Eigen::VectorXd& q, const Eigen::VectorXd& v, 
                             const HybridSolution& s, HybridDirection& d) const;

  void updateSolution(const std::vector<Robot>& robots, const HybridSolution& s,
                      const double primal_step_size, 
                      const double dual_step_size, HybridDirection& d);

  void computeKKTResidual(std::vector<Robot>& robots, 
                          const ContactSequence& contact_sequence,
                          const double t, const Eigen::VectorXd& q, 
                          const Eigen::VectorXd& v, const HybridSolution& s,
                          HybridKKTMatrix& kkt_matrix, 
                          HybridKKTResidual& kkt_residual);

  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

private:

  const Eigen::VectorXd& q_prev(const ContactSequence& contact_sequence, 
                                const HybridSolution& s,
                                const int time_stage) const;

  HybridOCP split_ocps_;
  TerminalOCP terminal_ocp_;
  double T_, dtau_;
  int N_, num_proc_;

};

} // namespace idocp 

#include "idocp/ocp/ocp_linearizer.hxx"

#endif // IDOCP_OCP_LINEARIZER_HPP_