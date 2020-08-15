#ifndef IDOCP_CONSTRAINTS_HPP_
#define IDOCP_CONSTRAINTS_HPP_

#include <vector>
#include <memory>

#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/constraints/constraint_component_base.hpp"
#include "idocp/constraints/constraints_data.hpp"


namespace idocp {

class Constraints {
public:
  Constraints() 
    : constraints_() {
  }

  ~Constraints() {}

  // Use default copy constructor.
  Constraints(const Constraints&) = default;

  // Use default copy coperator.
  Constraints& operator=(const Constraints&) = default;

  // Use default move constructor.
  Constraints(Constraints&&) noexcept = default;

  // Use default move assign coperator.
  Constraints& operator=(Constraints&&) noexcept = default;

  void push_back(const std::shared_ptr<ConstraintComponentBase>& constraint) {
    constraints_.push_back(constraint);
  }

  void clear() {
    constraints_.clear();
  }

  bool isEmpty() {
    return constraints_.empty();
  }

  ConstraintsData createConstraintsData(const Robot& robot) const {
    ConstraintsData datas;
    for (int i=0; i<constraints_.size(); ++i) {
      datas.data.push_back(ConstraintComponentData(constraints_[i]->dimc()));
    }
    return datas;
  }

  bool isFeasible(const Robot& robot, ConstraintsData& datas, 
                  const Eigen::Ref<const Eigen::VectorXd>& a, 
                  const Eigen::Ref<const Eigen::VectorXd>& f, 
                  const Eigen::Ref<const Eigen::VectorXd>& q, 
                  const Eigen::Ref<const Eigen::VectorXd>& v, 
                  const Eigen::Ref<const Eigen::VectorXd>& u) const {
    for (int i=0; i<constraints_.size(); ++i) {
      bool feasible = constraints_[i]->isFeasible(robot, datas.data[i], 
                                                  a, f, q, v, u);
      if (!feasible) {
        return false;
      }
    }
    return true;
  }

  void setSlackAndDual(const Robot& robot, ConstraintsData& datas, 
                       const double dtau, 
                       const Eigen::Ref<const Eigen::VectorXd>& a, 
                       const Eigen::Ref<const Eigen::VectorXd>& f, 
                       const Eigen::Ref<const Eigen::VectorXd>& q, 
                       const Eigen::Ref<const Eigen::VectorXd>& v, 
                       const Eigen::Ref<const Eigen::VectorXd>& u) const {
    for (int i=0; i<constraints_.size(); ++i) {
      constraints_[i]->setSlackAndDual(robot, datas.data[i], dtau, 
                                       a, f, q, v, u);
    }
  }

  void augmentDualResidual(const Robot& robot, ConstraintsData& datas,
                           const double dtau,
                           Eigen::Ref<Eigen::VectorXd> la, 
                           Eigen::Ref<Eigen::VectorXd> lf,
                           Eigen::Ref<Eigen::VectorXd> lq, 
                           Eigen::Ref<Eigen::VectorXd> lv) const {
    for (int i=0; i<constraints_.size(); ++i) {
      constraints_[i]->augmentDualResidual(robot, datas.data[i], dtau, 
                                           la, lf, lq, lv);
    }
  }

  void augmentDualResidual(const Robot& robot, ConstraintsData& datas,
                           const double dtau,
                           Eigen::Ref<Eigen::VectorXd> lu) const {
    for (int i=0; i<constraints_.size(); ++i) {
      constraints_[i]->augmentDualResidual(robot, datas.data[i], dtau, lu);
    }
  }


  void condenseSlackAndDual(const Robot& robot, ConstraintsData& datas,
                            const double dtau, 
                            const Eigen::Ref<const Eigen::VectorXd>& a, 
                            const Eigen::Ref<const Eigen::VectorXd>& f, 
                            const Eigen::Ref<const Eigen::VectorXd>& q, 
                            const Eigen::Ref<const Eigen::VectorXd>& v,
                            Eigen::Ref<Eigen::MatrixXd> Qaa, 
                            Eigen::Ref<Eigen::MatrixXd> Qff,
                            Eigen::Ref<Eigen::MatrixXd> Qqq, 
                            Eigen::Ref<Eigen::MatrixXd> Qvv,
                            Eigen::Ref<Eigen::VectorXd> la, 
                            Eigen::Ref<Eigen::VectorXd> lf,
                            Eigen::Ref<Eigen::VectorXd> lq, 
                            Eigen::Ref<Eigen::VectorXd> lv) const {
    for (int i=0; i<constraints_.size(); ++i) {
      constraints_[i]->condenseSlackAndDual(robot, datas.data[i], dtau, 
                                            a, f, q, v, Qaa, Qff, Qqq, Qvv,
                                            la, lf, lq, lv);
    }
  }

  void condenseSlackAndDual(const Robot& robot, ConstraintsData& datas,
                            const double dtau, 
                            const Eigen::Ref<const Eigen::VectorXd>& u, 
                            Eigen::Ref<Eigen::MatrixXd> Quu,
                            Eigen::Ref<Eigen::VectorXd> lu) const {
    for (int i=0; i<constraints_.size(); ++i) {
      constraints_[i]->condenseSlackAndDual(robot, datas.data[i], dtau, 
                                            u, Quu, lu);
    }
  }

  void computeSlackAndDualDirection(
      const Robot& robot, ConstraintsData& datas, const double dtau, 
      const Eigen::Ref<const Eigen::VectorXd>& da, 
      const Eigen::Ref<const Eigen::VectorXd>& df, 
      const Eigen::Ref<const Eigen::VectorXd>& dq, 
      const Eigen::Ref<const Eigen::VectorXd>& dv, 
      const Eigen::Ref<const Eigen::VectorXd>& du) const {
    for (int i=0; i<constraints_.size(); ++i) {
      constraints_[i]->computeSlackAndDualDirection(robot, datas.data[i], dtau, 
                                                    da, df, dq, dv, du);
    }
  }

  double maxSlackStepSize(const ConstraintsData& datas) const {
    double min_step_size = 1;
    for (int i=0; i<constraints_.size(); ++i) {
      const double step_size = constraints_[i]->maxSlackStepSize(datas.data[i]);
      if (step_size < min_step_size) {
        min_step_size = step_size;
      }
    }
    return min_step_size;
  }

  double maxDualStepSize(const ConstraintsData& datas) const {
    double min_step_size = 1;
    for (int i=0; i<constraints_.size(); ++i) {
      const double step_size = constraints_[i]->maxDualStepSize(datas.data[i]);
      if (step_size < min_step_size) {
        min_step_size = step_size;
      }
    }
    return min_step_size;
  }

  void updateSlack(ConstraintsData& datas, const double step_size) const {
    for (int i=0; i<constraints_.size(); ++i) {
      constraints_[i]->updateSlack(datas.data[i], step_size);
    }
  }

  void updateDual(ConstraintsData& datas, const double step_size) const {
    for (int i=0; i<constraints_.size(); ++i) {
      constraints_[i]->updateDual(datas.data[i], step_size);
    }
  }

  double costSlackBarrier(const ConstraintsData& datas) const {
    double cost = 0;
    for (int i=0; i<constraints_.size(); ++i) {
      cost += constraints_[i]->costSlackBarrier(datas.data[i]);
    }
    return cost;
  }

  double costSlackBarrier(const ConstraintsData& datas, 
                          const double step_size) const {
    double cost = 0;
    for (int i=0; i<constraints_.size(); ++i) {
      cost += constraints_[i]->costSlackBarrier(datas.data[i], step_size);
    }
    return cost;
  }

  double residualL1Nrom(const Robot& robot, ConstraintsData& datas, 
                        const double dtau, 
                        const Eigen::Ref<const Eigen::VectorXd>& a, 
                        const Eigen::Ref<const Eigen::VectorXd>& f, 
                        const Eigen::Ref<const Eigen::VectorXd>& q, 
                        const Eigen::Ref<const Eigen::VectorXd>& v, 
                        const Eigen::Ref<const Eigen::VectorXd>& u) const {
    double l1_norm = 0;
    for (int i=0; i<constraints_.size(); ++i) {
      l1_norm += constraints_[i]->residualL1Nrom(robot, datas.data[i], dtau,
                                                 a, f, q, v, u);
    }
    return l1_norm;
  }

  double residualSquaredNrom(const Robot& robot, ConstraintsData& datas, 
                             const double dtau, 
                             const Eigen::Ref<const Eigen::VectorXd>& a, 
                             const Eigen::Ref<const Eigen::VectorXd>& f, 
                             const Eigen::Ref<const Eigen::VectorXd>& q, 
                             const Eigen::Ref<const Eigen::VectorXd>& v, 
                             const Eigen::Ref<const Eigen::VectorXd>& u) const {
    double squared_norm = 0;
    for (int i=0; i<constraints_.size(); ++i) {
      squared_norm += constraints_[i]->residualSquaredNrom(robot, datas.data[i],
                                                           dtau, a, f, q, v, u);
    }
    return squared_norm;
  }


private:
  std::vector<std::shared_ptr<ConstraintComponentBase>> constraints_;

};

} // namespace idocp


#endif // IDOCP_CONSTRAINTS_HPP_