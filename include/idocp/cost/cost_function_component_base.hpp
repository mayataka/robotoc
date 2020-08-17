#ifndef IDOCP_COST_FUNCTION_COMPONENT_BASE_HPP_
#define IDOCP_COST_FUNCTION_COMPONENT_BASE_HPP_

#include "idocp/robot/robot.hpp"
#include "idocp/cost/cost_function_data.hpp"
#include "idocp/ocp/split_solution.hpp"
#include "idocp/ocp/kkt_residual.hpp"
#include "idocp/ocp/kkt_matrix.hpp"


namespace idocp {

class CostFunctionComponentBase {
public:
  CostFunctionComponentBase() {}

  virtual ~CostFunctionComponentBase() {}

  // Use default copy constructor.
  CostFunctionComponentBase(const CostFunctionComponentBase&) = default;

  // Use default copy coperator.
  CostFunctionComponentBase& operator=(const CostFunctionComponentBase&) 
      = default;

  // Use default move constructor.
  CostFunctionComponentBase(CostFunctionComponentBase&&) noexcept = default;

  // Use default move assign coperator.
  CostFunctionComponentBase& operator=(CostFunctionComponentBase&&) noexcept 
      = default;

  virtual double l(const Robot& robot, CostFunctionData& data, const double t, 
                   const double dtau, const SplitSolution& s) const = 0;

  virtual double phi(const Robot& robot, CostFunctionData& data, const double t, 
                     const SplitSolution& s) const = 0;

  virtual void lq(const Robot& robot, CostFunctionData& data, const double t, 
                  const double dtau, const SplitSolution& s, 
                  KKTResidual& kkt_residual) const = 0;

  virtual void lv(const Robot& robot, CostFunctionData& data, const double t, 
                  const double dtau, const SplitSolution& s, 
                  KKTResidual& kkt_residual) const = 0;

  virtual void la(const Robot& robot, CostFunctionData& data, const double t, 
                  const double dtau, const SplitSolution& s, 
                  KKTResidual& kkt_residual) const = 0;

  virtual void lf(const Robot& robot, CostFunctionData& data, const double t, 
                  const double dtau, const SplitSolution& s, 
                  KKTResidual& kkt_residual) const = 0;

  virtual void lu(const Robot& robot, CostFunctionData& data, const double t, 
                  const double dtau, const SplitSolution& s, 
                  KKTResidual& kkt_residual) const = 0;

  virtual void lqq(const Robot& robot, CostFunctionData& data, const double t, 
                   const double dtau, const SplitSolution& s, 
                   KKTMatrix& kkt_matrix) const = 0;

  virtual void lvv(const Robot& robot, CostFunctionData& data, const double t, 
                   const double dtau, const SplitSolution& s, 
                   KKTMatrix& kkt_matrix) const = 0;

  virtual void laa(const Robot& robot, CostFunctionData& data, const double t, 
                   const double dtau, const SplitSolution& s, 
                   KKTMatrix& kkt_matrix) const = 0;

  virtual void lff(const Robot& robot, CostFunctionData& data, const double t, 
                   const double dtau, const SplitSolution& s, 
                   KKTMatrix& kkt_matrix) const = 0;

  virtual void luu(const Robot& robot, CostFunctionData& data, const double t, 
                   const double dtau, const SplitSolution& s, 
                   KKTMatrix& kkt_matrix) const = 0;

  virtual void phiq(const Robot& robot, CostFunctionData& data,  
                    const double t, const SplitSolution& s, 
                    KKTResidual& kkt_residual) const = 0;

  virtual void phiv(const Robot& robot, CostFunctionData& data,  
                    const double t, const SplitSolution& s, 
                    KKTResidual& kkt_residual) const = 0;

  virtual void phiqq(const Robot& robot, CostFunctionData& data,  
                     const double t, const SplitSolution& s,
                     KKTMatrix& kkt_matrix) const = 0;

  virtual void phivv(const Robot& robot, CostFunctionData& data,  
                     const double t, const SplitSolution& s,
                     KKTMatrix& kkt_matrix) const = 0;

};

} // namespace idocp

#endif // IDOCP_COST_FUNCTION_COMPONENT_BASE_HPP_