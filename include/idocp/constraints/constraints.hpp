#ifndef IDOCP_CONSTRAINTS_HPP_
#define IDOCP_CONSTRAINTS_HPP_

#include <vector>
#include <memory>

#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/ocp/split_solution.hpp"
#include "idocp/ocp/split_direction.hpp"
#include "idocp/constraints/constraint_component_base.hpp"
#include "idocp/constraints/constraint_component_data.hpp"
#include "idocp/ocp/kkt_residual.hpp"
#include "idocp/ocp/kkt_matrix.hpp"


namespace idocp {

///
/// @typedef ConstraintsData
/// @brief Data for constraints. Composed of ConstraintComponentData 
/// corrensponding to the components of Constraints.
///
typedef std::vector<ConstraintComponentData> ConstraintsData;

class Constraints {
public:
  Constraints();

  ~Constraints();

  // Use default copy constructor.
  Constraints(const Constraints&) = default;

  // Use default copy coperator.
  Constraints& operator=(const Constraints&) = default;

  // Use default move constructor.
  Constraints(Constraints&&) noexcept = default;

  // Use default move assign coperator.
  Constraints& operator=(Constraints&&) noexcept = default;

  void push_back(const std::shared_ptr<ConstraintComponentBase>& constraint);

  void clear();

  bool isEmpty() const;

  bool useKinematics() const;

  ConstraintsData createConstraintsData(const Robot& robot) const;

  bool isFeasible(Robot& robot, ConstraintsData& data,
                  const SplitSolution& s) const;

  void setSlackAndDual(Robot& robot, ConstraintsData& data, 
                       const double dtau, const SplitSolution& s) const;

  void augmentDualResidual(Robot& robot, ConstraintsData& data,
                           const double dtau, KKTResidual& kkt_residual) const;

  void augmentDualResidual(const Robot& robot, ConstraintsData& data,
                           const double dtau, Eigen::VectorXd& lu) const;

  void condenseSlackAndDual(Robot& robot, ConstraintsData& data,
                            const double dtau, const SplitSolution& s,
                            KKTMatrix& kkt_matrix, 
                            KKTResidual& kkt_residual) const;

  void condenseSlackAndDual(const Robot& robot, ConstraintsData& data,
                            const double dtau, const Eigen::VectorXd& u,
                            Eigen::MatrixXd& Quu, Eigen::VectorXd& lu) const;

  void computeSlackAndDualDirection(Robot& robot, ConstraintsData& data, 
                                    const double dtau, 
                                    const SplitDirection& d) const;

  double maxSlackStepSize(const ConstraintsData& data) const;

  double maxDualStepSize(const ConstraintsData& data) const;

  void updateSlack(ConstraintsData& data, const double step_size) const;

  void updateDual(ConstraintsData& data, const double step_size) const;

  double costSlackBarrier(const ConstraintsData& data) const;

  double costSlackBarrier(const ConstraintsData& data, 
                          const double step_size) const;

  double residualL1Nrom(Robot& robot, ConstraintsData& data, 
                        const double dtau, const SplitSolution& s) const;

  double squaredKKTErrorNorm(Robot& robot, ConstraintsData& data, 
                             const double dtau, const SplitSolution& s) const;

private:
  std::vector<std::shared_ptr<ConstraintComponentBase>> constraints_;

};

} // namespace idocp

#include "idocp/constraints/constraints.hxx"

#endif // IDOCP_CONSTRAINTS_HPP_