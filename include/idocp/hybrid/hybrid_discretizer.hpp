#ifndef IDOCP_HYBRID_DISCRETIZER_HPP_
#define IDOCP_HYBRID_DISCRETIZER_HPP_

#include "idocp/robot/robot.hpp"

#include "idocp/ocp/split_ocp.hpp"
#include "idocp/impulse/split_impulse_ocp.hpp"
#include "idocp/ocp/terminal_ocp.hpp"
#include "idocp/cost/cost_function.hpp"
#include "idocp/constraints/constraints.hpp"
#include "idocp/ocp/split_solution.hpp"
#include "idocp/ocp/split_direction.hpp"
#include "idocp/ocp/split_kkt_matrix.hpp"
#include "idocp/ocp/split_kkt_residual.hpp"
#include "idocp/impulse/impulse_split_solution.hpp"
#include "idocp/impulse/impulse_split_direction.hpp"
#include "idocp/impulse/impulse_split_kkt_matrix.hpp"
#include "idocp/impulse/impulse_split_kkt_residual.hpp"
#include "idocp/ocp/riccati_factorization.hpp"
#include "idocp/ocp/riccati_factorizer.hpp"
#include "idocp/impulse/impulse_riccati_factorizer.hpp"

#include <vector>
#include <memory>
#include <cassert>

namespace idocp {

class HybridDiscretizer {
public:
  HybridDiscretizer();

  int N_total() const;

  int N_impulse() const;

  int N_lift() const;

  int

  double dtau(const int time_stage) const;

  double dtau_after_impulse(const int impulse_index) const;

  double dtau_after_lift(const int lift_index) const;

  const Eigen::VectorXd& q_prev() const;

};

} // namespace idocp

#endif // IDOCP_HYBRID_DISCRETIZER_HPP_ 