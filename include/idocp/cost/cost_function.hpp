#ifndef IDOCP_COST_FUNCTION_HPP_
#define IDOCP_COST_FUNCTION_HPP_

#include <vector>
#include <memory>
#include <assert.h>

#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/cost/cost_function_component_base.hpp"
#include "idocp/cost/cost_function_data.hpp"


namespace idocp {

class CostFunction {
public:
  CostFunction() 
    : costs_() {
  }

  ~CostFunction() {}

  // Use default copy constructor.
  CostFunction(const CostFunction&) = default;

  // Use default copy coperator.
  CostFunction& operator=(const CostFunction&) = default;

  // Use default move constructor.
  CostFunction(CostFunction&&) noexcept = default;

  // Use default move assign coperator.
  CostFunction& operator=(CostFunction&&) noexcept = default;

  void push_back(const std::shared_ptr<CostFunctionComponentBase>& cost) {
    costs_.push_back(cost);
  }

  void clear() {
    costs_.clear();
  }

  bool isEmpty() {
    return costs_.empty();
  }

  double l(const Robot& robot, CostFunctionData& data, const double t, 
           const double dtau, const Eigen::Ref<const Eigen::VectorXd>& q, 
           const Eigen::Ref<const Eigen::VectorXd>& v, 
           const Eigen::Ref<const Eigen::VectorXd>& a, 
           const Eigen::Ref<const Eigen::VectorXd>& f, 
           const Eigen::Ref<const Eigen::VectorXd>& u) const {
    assert(dtau > 0);
    assert(q.size() == robot.dimq());
    assert(v.size() == robot.dimv());
    assert(a.size() == robot.dimv());
    assert(f.size() >= robot.dimf());
    assert(u.size() == robot.dimv());
    double l = 0;
    for (const auto cost : costs_) {
      l += cost->l(robot, data, t, dtau, q, v, a, f, u);
    }
    return l;
  }

  double phi(const Robot& robot, CostFunctionData& data, const double t, 
             const Eigen::Ref<const Eigen::VectorXd>& q, 
             const Eigen::Ref<const Eigen::VectorXd>& v) const {
    assert(q.size() == robot.dimq());
    assert(v.size() == robot.dimv());
    double phi = 0;
    for (const auto cost : costs_) {
      phi += cost->phi(robot, data, t, q, v);
    }
    return phi;
  }

  void lq(const Robot& robot, CostFunctionData& data, const double t, 
          const double dtau, const Eigen::Ref<const Eigen::VectorXd>& q, 
          const Eigen::Ref<const Eigen::VectorXd>& v, 
          const Eigen::Ref<const Eigen::VectorXd>& a, 
          Eigen::Ref<Eigen::VectorXd> lq) const {
    assert(dtau > 0);
    assert(q.size() == robot.dimq());
    assert(v.size() == robot.dimv());
    assert(a.size() == robot.dimv());
    assert(lq.size() == robot.dimv());
    for (const auto cost : costs_) {
      cost->lq(robot, data, t, dtau, q, v, a, lq);
    }
  }

  void lv(const Robot& robot, CostFunctionData& data, const double t, 
          const double dtau, const Eigen::Ref<const Eigen::VectorXd>& q, 
          const Eigen::Ref<const Eigen::VectorXd>& v, 
          const Eigen::Ref<const Eigen::VectorXd>& a, 
          Eigen::Ref<Eigen::VectorXd> lv) const {
    assert(dtau > 0);
    assert(q.size() == robot.dimq());
    assert(v.size() == robot.dimv());
    assert(a.size() == robot.dimv());
    assert(lv.size() == robot.dimv());
    for (const auto cost : costs_) {
      cost->lv(robot, data, t, dtau, q, v, a, lv);
    }
  }

  void la(const Robot& robot, CostFunctionData& data, const double t, 
          const double dtau, const Eigen::Ref<const Eigen::VectorXd>& q, 
          const Eigen::Ref<const Eigen::VectorXd>& v, 
          const Eigen::Ref<const Eigen::VectorXd>& a, 
          Eigen::Ref<Eigen::VectorXd> la) const {
    assert(dtau > 0);
    assert(q.size() == robot.dimq());
    assert(v.size() == robot.dimv());
    assert(a.size() == robot.dimv());
    assert(la.size() == robot.dimv());
    for (const auto cost : costs_) {
      cost->la(robot, data, t, dtau, q, v, a, la);
    }
  }

  void lf(const Robot& robot, CostFunctionData& data, const double t, 
          const double dtau, const Eigen::Ref<const Eigen::VectorXd>& f, 
          Eigen::Ref<Eigen::VectorXd> lf) const {
    assert(dtau > 0);
    assert(f.size() >= robot.dimf());
    assert(lf.size() >= robot.dimf());
    for (const auto cost : costs_) {
      cost->lf(robot, data, t, dtau, f, lf);
    }
  }

  void lu(const Robot& robot, CostFunctionData& data, const double t, 
          const double dtau, const Eigen::Ref<const Eigen::VectorXd>& u, 
          Eigen::Ref<Eigen::VectorXd> lu) const {
    assert(dtau > 0);
    assert(u.size() == robot.dimv());
    assert(lu.size() == robot.dimv());
    for (const auto cost : costs_) {
      cost->lu(robot, data, t, dtau, u, lu);
    }
  }

  void lqq(const Robot& robot, CostFunctionData& data, const double t, 
           const double dtau, const Eigen::Ref<const Eigen::VectorXd>& q, 
           const Eigen::Ref<const Eigen::VectorXd>& v, 
           const Eigen::Ref<const Eigen::VectorXd>& a, 
           Eigen::Ref<Eigen::MatrixXd> lqq) const {
    assert(dtau > 0);
    assert(q.size() == robot.dimq());
    assert(v.size() == robot.dimv());
    assert(a.size() == robot.dimv());
    assert(lqq.rows() == robot.dimv());
    assert(lqq.cols() == robot.dimv());
    for (const auto cost : costs_) {
      cost->lqq(robot, data, t, dtau, q, v, a, lqq);
    }
  }

  void lvv(const Robot& robot, CostFunctionData& data, const double t, 
           const double dtau, const Eigen::Ref<const Eigen::VectorXd>& q, 
           const Eigen::Ref<const Eigen::VectorXd>& v, 
           const Eigen::Ref<const Eigen::VectorXd>& a, 
           Eigen::Ref<Eigen::MatrixXd> lvv) const {
    assert(dtau > 0);
    assert(q.size() == robot.dimq());
    assert(v.size() == robot.dimv());
    assert(a.size() == robot.dimv());
    assert(lvv.rows() == robot.dimv());
    assert(lvv.cols() == robot.dimv());
    for (const auto cost : costs_) {
      cost->lvv(robot, data, t, dtau, q, v, a, lvv);
    }
  }

  void laa(const Robot& robot, CostFunctionData& data, const double t, 
           const double dtau, const Eigen::Ref<const Eigen::VectorXd>& q, 
           const Eigen::Ref<const Eigen::VectorXd>& v, 
           const Eigen::Ref<const Eigen::VectorXd>& a, 
           Eigen::Ref<Eigen::MatrixXd> laa) const {
    assert(dtau > 0);
    assert(q.size() == robot.dimq());
    assert(v.size() == robot.dimv());
    assert(a.size() == robot.dimv());
    assert(laa.rows() == robot.dimv());
    assert(laa.cols() == robot.dimv());
    for (const auto cost : costs_) {
      cost->laa(robot, data, t, dtau, q, v, a, laa);
    }
  }

  void lff(const Robot& robot, CostFunctionData& data, const double t, 
           const double dtau, const Eigen::Ref<const Eigen::VectorXd>& f, 
           Eigen::Ref<Eigen::MatrixXd> lff) const {
    assert(dtau > 0);
    assert(f.size() >= robot.dimf());
    assert(lff.rows() >= robot.dimf());
    assert(lff.cols() >= robot.dimf());
    for (const auto cost : costs_) {
      cost->lff(robot, data, t, dtau, f, lff);
    }
  }

  void luu(const Robot& robot, CostFunctionData& data, const double t, 
           const double dtau, const Eigen::Ref<const Eigen::VectorXd>& u, 
           Eigen::Ref<Eigen::MatrixXd> luu) const {
    assert(dtau > 0);
    assert(u.size() == robot.dimv());
    assert(luu.rows() == robot.dimv());
    assert(luu.cols() == robot.dimv());
    for (const auto cost : costs_) {
      cost->luu(robot, data, t, dtau, u, luu);
    }
  }

  void augment_lqq(const Robot& robot, CostFunctionData& data, 
                   const double t, const double dtau, 
                   const Eigen::Ref<const Eigen::VectorXd>& q, 
                   const Eigen::Ref<const Eigen::VectorXd>& v, 
                   const Eigen::Ref<const Eigen::VectorXd>& a, 
                   Eigen::Ref<Eigen::MatrixXd> lqq) const {
    assert(dtau > 0);
    assert(q.size() == robot.dimq());
    assert(v.size() == robot.dimv());
    assert(a.size() == robot.dimv());
    assert(lqq.rows() == robot.dimv());
    assert(lqq.cols() == robot.dimv());
    for (const auto cost : costs_) {
      cost->augment_lqq(robot, data, t, dtau, q, v, a, lqq);
    }
  }

  void augment_lvv(const Robot& robot, CostFunctionData& data, 
                   const double t, const double dtau, 
                   const Eigen::Ref<const Eigen::VectorXd>& q, 
                   const Eigen::Ref<const Eigen::VectorXd>& v, 
                   const Eigen::Ref<const Eigen::VectorXd>& a, 
                   Eigen::Ref<Eigen::MatrixXd> lvv) const {
    assert(dtau > 0);
    assert(q.size() == robot.dimq());
    assert(v.size() == robot.dimv());
    assert(a.size() == robot.dimv());
    assert(lvv.rows() == robot.dimv());
    assert(lvv.cols() == robot.dimv());
    for (const auto cost : costs_) {
      cost->augment_lvv(robot, data, t, dtau, q, v, a, lvv);
    }
  }

  void augment_laa(const Robot& robot, CostFunctionData& data, 
                   const double t, const double dtau, 
                   const Eigen::Ref<const Eigen::VectorXd>& q, 
                   const Eigen::Ref<const Eigen::VectorXd>& v, 
                   const Eigen::Ref<const Eigen::VectorXd>& a, 
                   Eigen::Ref<Eigen::MatrixXd> laa) const {
    assert(dtau > 0);
    assert(q.size() == robot.dimq());
    assert(v.size() == robot.dimv());
    assert(a.size() == robot.dimv());
    assert(laa.rows() == robot.dimv());
    assert(laa.cols() == robot.dimv());
    for (const auto cost : costs_) {
      cost->augment_laa(robot, data, t, dtau, q, v, a, laa);
    }
  }

  void augment_lff(const Robot& robot, CostFunctionData& data, 
                   const double t, const double dtau, 
                   const Eigen::Ref<const Eigen::VectorXd>& f, 
                   Eigen::Ref<Eigen::MatrixXd> lff) const {
    assert(dtau > 0);
    assert(f.size() >= robot.dimf());
    assert(lff.rows() >= robot.dimf());
    assert(lff.cols() >= robot.dimf());
    for (const auto cost : costs_) {
      cost->augment_lff(robot, data, t, dtau, f, lff);
    }
  }

  void augment_luu(const Robot& robot, CostFunctionData& data, 
                   const double t, const double dtau, 
                   const Eigen::Ref<const Eigen::VectorXd>& u, 
                   Eigen::Ref<Eigen::MatrixXd> luu) const {
    assert(dtau > 0);
    assert(u.size() == robot.dimv());
    assert(luu.rows() == robot.dimv());
    assert(luu.cols() == robot.dimv());
    for (const auto cost : costs_) {
      cost->augment_luu(robot, data, t, dtau, u, luu);
    }
  }

  void phiq(const Robot& robot, CostFunctionData& data, const double t, 
            const Eigen::Ref<const Eigen::VectorXd>& q, 
            const Eigen::Ref<const Eigen::VectorXd>& v, 
            Eigen::Ref<Eigen::VectorXd> phiq) const {
    assert(q.size() == robot.dimq());
    assert(v.size() == robot.dimv());
    assert(phiq.size() == robot.dimv());
    for (const auto cost : costs_) {
      cost->phiq(robot, data, t, q, v, phiq);
    }
  }

  void phiv(const Robot& robot, CostFunctionData& data, const double t, 
            const Eigen::Ref<const Eigen::VectorXd>& q, 
            const Eigen::Ref<const Eigen::VectorXd>& v, 
            Eigen::Ref<Eigen::VectorXd> phiv) const {
    assert(q.size() == robot.dimq());
    assert(v.size() == robot.dimv());
    assert(phiv.size() == robot.dimv());
    for (const auto cost : costs_) {
      cost->phiv(robot, data, t, q, v, phiv);
    }
  }

  void phiqq(const Robot& robot, CostFunctionData& data, const double t, 
             const Eigen::Ref<const Eigen::VectorXd>& q, 
             const Eigen::Ref<const Eigen::VectorXd>& v, 
             Eigen::Ref<Eigen::MatrixXd> phiqq) const {
    assert(q.size() == robot.dimq());
    assert(v.size() == robot.dimv());
    assert(phiqq.rows() == robot.dimv());
    assert(phiqq.cols() == robot.dimv());
    for (const auto cost : costs_) {
      cost->phiqq(robot, data, t, q, v, phiqq);
    }
  }

  void phivv(const Robot& robot, CostFunctionData& data, const double t, 
             const Eigen::Ref<const Eigen::VectorXd>& q, 
             const Eigen::Ref<const Eigen::VectorXd>& v, 
             Eigen::Ref<Eigen::MatrixXd> phivv) const {
    assert(q.size() == robot.dimq());
    assert(v.size() == robot.dimv());
    assert(phivv.rows() == robot.dimv());
    assert(phivv.cols() == robot.dimv());
    for (const auto cost : costs_) {
      cost->phivv(robot, data, t, q, v, phivv);
    }
  }

  void augment_phiqq(const Robot& robot, CostFunctionData& data, const double t, 
                     const Eigen::Ref<const Eigen::VectorXd>& q, 
                     const Eigen::Ref<const Eigen::VectorXd>& v, 
                     Eigen::Ref<Eigen::MatrixXd> phiqq) const {
    assert(q.size() == robot.dimq());
    assert(v.size() == robot.dimv());
    assert(phiqq.rows() == robot.dimv());
    assert(phiqq.cols() == robot.dimv());
    for (const auto cost : costs_) {
      cost->augment_phiqq(robot, data, t, q, v, phiqq);
    }
  }

  void augment_phivv(const Robot& robot, CostFunctionData& data, const double t, 
                     const Eigen::Ref<const Eigen::VectorXd>& q, 
                     const Eigen::Ref<const Eigen::VectorXd>& v, 
                     Eigen::Ref<Eigen::MatrixXd> phivv) const {
    assert(q.size() == robot.dimq());
    assert(v.size() == robot.dimv());
    assert(phivv.rows() == robot.dimv());
    assert(phivv.cols() == robot.dimv());
    for (const auto cost : costs_) {
      cost->augment_phivv(robot, data, t, q, v, phivv);
    }
  }

private:
  std::vector<std::shared_ptr<CostFunctionComponentBase>> costs_;

};

} // namespace idocp


#endif // IDOCP_COST_FUNCTION_HPP_