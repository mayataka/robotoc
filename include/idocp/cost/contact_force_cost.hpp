#ifndef IDOCP_CONTACT_FORCE_COST_HPP_
#define IDOCP_CONTACT_FORCE_COST_HPP_

#include <vector>

#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/cost/cost_function_component_base.hpp"
#include "idocp/cost/cost_function_data.hpp"
#include "idocp/ocp/split_solution.hpp"
#include "idocp/ocp/split_kkt_residual.hpp"
#include "idocp/ocp/split_kkt_matrix.hpp"


namespace idocp {

///
/// @class ContactForceCost
/// @brief Cost on the contact forces expressed in the local frames.
///
class ContactForceCost final : public CostFunctionComponentBase {
public:
  ///
  /// @brief Constructor. 
  /// @param[in] robot Robot model.
  ///
  ContactForceCost(const Robot& robot);

  ///
  /// @brief Default constructor. 
  ///
  ContactForceCost();

  ///
  /// @brief Destructor. 
  ///
  ~ContactForceCost();

  ///
  /// @brief Default copy constructor. 
  ///
  ContactForceCost(const ContactForceCost&) = default;

  ///
  /// @brief Default copy operator. 
  ///
  ContactForceCost& operator=(const ContactForceCost&) = default;

  ///
  /// @brief Default move constructor. 
  ///
  ContactForceCost(ContactForceCost&&) noexcept = default;

  ///
  /// @brief Default move assign operator. 
  ///
  ContactForceCost& operator=(ContactForceCost&&) noexcept = default;

  ///
  /// @brief Sets the reference contact forces expressed in the local frames. 
  /// @param[in] f_ref Reference contact forces expressed in the local frames. 
  /// Size must be Robot::maxPointContacts().
  ///
  void set_f_ref(const std::vector<Eigen::Vector3d>& f_ref);

  ///
  /// @brief Sets the weight vectors on the contact forces. 
  /// @param[in] f_weight Weight vectors on the contact forces. 
  /// Size must be Robot::maxPointContacts().
  ///
  void set_f_weight(const std::vector<Eigen::Vector3d>& f_weight);

  ///
  /// @brief Sets the reference impulse forces expressed in the local frames. 
  /// @param[in] fi_ref Reference impulse forces expressed in the local frames. 
  /// Size must be Robot::maxPointContacts().
  ///
  void set_fi_ref(const std::vector<Eigen::Vector3d>& fi_ref);

  ///
  /// @brief Sets the weight vectors on the impulse forces. 
  /// @param[in] fi_weight Weight vectors on the impulse forces. 
  /// Size must be Robot::maxPointContacts().
  ///
  void set_fi_weight(const std::vector<Eigen::Vector3d>& fi_weight);

  bool useKinematics() const override;

  double computeStageCost(Robot& robot, CostFunctionData& data, const double t, 
                          const double dt, 
                          const SplitSolution& s) const override;

  void computeStageCostDerivatives(
      Robot& robot, CostFunctionData& data, const double t, const double dt, 
      const SplitSolution& s, SplitKKTResidual& kkt_residual) const override;

  void computeStageCostHessian(Robot& robot, CostFunctionData& data, 
                               const double t, const double dt, 
                               const SplitSolution& s, 
                               SplitKKTMatrix& kkt_matrix) const override;

  double computeTerminalCost(Robot& robot, CostFunctionData& data, 
                             const double t, 
                             const SplitSolution& s) const override;

  void computeTerminalCostDerivatives(
      Robot& robot, CostFunctionData& data, const double t, 
      const SplitSolution& s, SplitKKTResidual& kkt_residual) const override;

  void computeTerminalCostHessian(Robot& robot, CostFunctionData& data, 
                                  const double t, const SplitSolution& s, 
                                  SplitKKTMatrix& kkt_matrix) const override;

  double computeImpulseCost(Robot& robot, CostFunctionData& data, 
                            const double t, 
                            const ImpulseSplitSolution& s) const override;

  void computeImpulseCostDerivatives(
      Robot& robot, CostFunctionData& data, const double t, 
      const ImpulseSplitSolution& s, 
      ImpulseSplitKKTResidual& kkt_residual) const;

  void computeImpulseCostHessian(
      Robot& robot, CostFunctionData& data, const double t, 
      const ImpulseSplitSolution& s, 
      ImpulseSplitKKTMatrix& kkt_matrix) const override;

  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

private:
  int max_point_contacts_, max_dimf_;
  std::vector<Eigen::Vector3d> f_ref_, f_weight_, fi_ref_, fi_weight_;

};

} // namespace idocp


#endif // IDOCP_CONTACT_FORCE_COST_HPP_ 