#include "idocp/constraints/constraints.hpp"


namespace idocp {

Constraints::Constraints() 
  : constraints_() {
}


Constraints::~Constraints() {
}


void Constraints::push_back(
    const std::shared_ptr<ConstraintComponentBase>& constraint) {
  constraints_.push_back(constraint);
}


void Constraints::clear() {
  constraints_.clear();
}


bool Constraints::isEmpty() const {
  return constraints_.empty();
}


ConstraintsData Constraints::createConstraintsData(const Robot& robot) const {
  ConstraintsData datas;
  for (int i=0; i<constraints_.size(); ++i) {
    datas.data.push_back(ConstraintComponentData(constraints_[i]->dimc()));
  }
  return datas;
}


bool Constraints::isFeasible(const Robot& robot, ConstraintsData& datas, 
                             const SplitSolution& s) const {
  for (int i=0; i<constraints_.size(); ++i) {
    bool feasible = constraints_[i]->isFeasible(robot, datas.data[i], s);
    if (!feasible) {
      return false;
    }
  }
  return true;
}


void Constraints::setSlackAndDual(const Robot& robot, ConstraintsData& datas, 
                                  const double dtau, 
                                  const SplitSolution& s) const {
  for (int i=0; i<constraints_.size(); ++i) {
    constraints_[i]->setSlackAndDual(robot, datas.data[i], dtau, s);
  }
}


void Constraints::augmentDualResidual(const Robot& robot, 
                                      ConstraintsData& datas, const double dtau, 
                                      KKTResidual& kkt_residual) const {
  for (int i=0; i<constraints_.size(); ++i) {
    constraints_[i]->augmentDualResidual(robot, datas.data[i], dtau, 
                                          kkt_residual);
  }
}


void Constraints::condenseSlackAndDual(const Robot& robot, 
                                       ConstraintsData& datas, 
                                       const double dtau, 
                                       const SplitSolution& s, 
                                       KKTMatrix& kkt_matrix, 
                                       KKTResidual& kkt_residual) const {
  for (int i=0; i<constraints_.size(); ++i) {
    constraints_[i]->condenseSlackAndDual(robot, datas.data[i], dtau, s,
                                          kkt_matrix, kkt_residual);
  }
}


void Constraints::computeSlackAndDualDirection(const Robot& robot, 
                                               ConstraintsData& datas, 
                                               const double dtau, 
                                               const SplitDirection& d) const {
  for (int i=0; i<constraints_.size(); ++i) {
    constraints_[i]->computeSlackAndDualDirection(robot, datas.data[i], dtau, d);
  }
}


double Constraints::maxSlackStepSize(const ConstraintsData& datas) const {
  double min_step_size = 1;
  for (int i=0; i<constraints_.size(); ++i) {
    const double step_size = constraints_[i]->maxSlackStepSize(datas.data[i]);
    if (step_size < min_step_size) {
      min_step_size = step_size;
    }
  }
  return min_step_size;
}


double Constraints::maxDualStepSize(const ConstraintsData& datas) const {
  double min_step_size = 1;
  for (int i=0; i<constraints_.size(); ++i) {
    const double step_size = constraints_[i]->maxDualStepSize(datas.data[i]);
    if (step_size < min_step_size) {
      min_step_size = step_size;
    }
  }
  return min_step_size;
}


void Constraints::updateSlack(ConstraintsData& datas, 
                              const double step_size) const {
  for (int i=0; i<constraints_.size(); ++i) {
    constraints_[i]->updateSlack(datas.data[i], step_size);
  }
}


void Constraints::updateDual(ConstraintsData& datas, 
                             const double step_size) const {
  for (int i=0; i<constraints_.size(); ++i) {
    constraints_[i]->updateDual(datas.data[i], step_size);
  }
}


double Constraints::costSlackBarrier(const ConstraintsData& datas) const {
  double cost = 0;
  for (int i=0; i<constraints_.size(); ++i) {
    cost += constraints_[i]->costSlackBarrier(datas.data[i]);
  }
  return cost;
}


double Constraints::costSlackBarrier(const ConstraintsData& datas, 
                        const double step_size) const {
  double cost = 0;
  for (int i=0; i<constraints_.size(); ++i) {
    cost += constraints_[i]->costSlackBarrier(datas.data[i], step_size);
  }
  return cost;
}


double Constraints::residualL1Nrom(const Robot& robot, ConstraintsData& datas, 
                                   const double dtau, 
                                   const SplitSolution& s) const {
  double l1_norm = 0;
  for (int i=0; i<constraints_.size(); ++i) {
    l1_norm += constraints_[i]->residualL1Nrom(robot, datas.data[i], dtau, s);
  }
  return l1_norm;
}


double Constraints::squaredKKTErrorNorm(const Robot& robot, 
                                        ConstraintsData& datas, 
                                        const double dtau, 
                                        const SplitSolution& s) const {
  double squared_norm = 0;
  for (int i=0; i<constraints_.size(); ++i) {
    squared_norm += constraints_[i]->squaredKKTErrorNorm(robot, datas.data[i],
                                                          dtau, s);
  }
  return squared_norm;
}

} // namespace idocp