#ifndef IDImpulseOCP_IMPULSE_ImpulseOCP_HPP_
#define IDImpulseOCP_IMPULSE_ImpulseOCP_HPP_

#include <vector>
#include <memory>

#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/cost/cost_function.hpp"
#include "idocp/constraints/constraints.hpp"
#include "idocp/cost/impulse_cost_function.hpp"
#include "idocp/constraints/impulse_constraints.hpp"
#include "idocp/ocp/contact_sequence.hpp"
#include "idocp/ocp/split_ocp.hpp"
#include "idocp/impulse/split_impulse_ocp.hpp"
#include "idocp/ocp/terminal_ocp.hpp"
#include "idocp/ocp/split_solution.hpp"
#include "idocp/ocp/split_direction.hpp"
#include "idocp/impulse/impulse_split_solution.hpp"
#include "idocp/impulse/impulse_split_direction.hpp"
#include "idocp/ocp/riccati_solution.hpp"
#include "idocp/ocp/line_search_filter.hpp"
#include "idocp/ocp/split_temporary_solution.hpp"


namespace idocp {

///
/// @class ImpulseOCP
/// @brief Optimal control problem solver by Riccati recursion. 
///
class ImpulseOCP {
public:
  ///
  /// @brief Construct optimal control problem solver.
  /// @param[in] robot Robot model. Must be initialized by URDF or XML.
  /// @param[in] cost Shared ptr to the cost function.
  /// @param[in] constraints Shared ptr to the constraints.
  ///
  ImpulseOCP(const Robot& robot, 
             const std::shared_ptr<ImpulseCostFunction>& cost, 
             const std::shared_ptr<ImpulseConstraints>& constraints);

  ///
  /// @brief Default constructor. 
  ///
  ImpulseOCP();

  ///
  /// @brief Destructor. 
  ///
  ~ImpulseOCP();

  ///
  /// @brief Default copy constructor. 
  ///
  ImpulseOCP(const ImpulseOCP&) = default;

  ///
  /// @brief Default copy assign operator. 
  ///
  ImpulseOCP& operator=(const ImpulseOCP&) = default;

  ///
  /// @brief Default move constructor. 
  ///
  ImpulseOCP(ImpulseOCP&&) noexcept = default;

  ///
  /// @brief Default move assign operator. 
  ///
  ImpulseOCP& operator=(ImpulseOCP&&) noexcept = default;

  void formOCP(const ContactSequence& contact_sequence);

  SplitImpulseOCP& splitOCP(const int time_stage);

  const SplitImpulseOCP& splitOCP(const int time_stage) const;

  ImpulseSplitSolution& s(const int stage);

  const ImpulseSplitSolution& s(const int stage) const;

  ImpulseSplitDirection& d(const int stage);

  const ImpulseSplitDirection& d(const int stage) const;

  RiccatiFactorization& riccati(const int stage);

  const RiccatiFactorization& riccati(const int stage) const;

  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

private:
  std::vector<SplitImpulseImpulseOCP> split_impulse_ocps_; 
  std::vector<ImpulseSplitSolution> s_; 
  std::vector<ImpulseSplitDirection> d_; 
  std::vector<RiccatiFactorization> riccati_;
};

} // namespace idocp 


#endif // IDImpulseOCP_IMPULSE_ImpulseOCP_HPP_ 