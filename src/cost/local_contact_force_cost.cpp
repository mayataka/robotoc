#include "robotoc/cost/local_contact_force_cost.hpp"

#include <iostream>
#include <stdexcept>


namespace robotoc {

LocalContactForceCost::LocalContactForceCost(const Robot& robot)
  : CostFunctionComponentBase(),
    max_num_contacts_(robot.maxNumContacts()),
    max_dimf_(robot.max_dimf()),
    contact_types_(robot.contactTypes()),
    f_ref_(robot.maxNumContacts(), Eigen::Vector3d::Zero()),
    f_weight_(robot.maxNumContacts(), Eigen::Vector3d::Zero()),
    fi_ref_(robot.maxNumContacts(), Eigen::Vector3d::Zero()),
    fi_weight_(robot.maxNumContacts(), Eigen::Vector3d::Zero()) {
}


LocalContactForceCost::LocalContactForceCost()
  : CostFunctionComponentBase(),
    max_num_contacts_(0),
    max_dimf_(0),
    contact_types_(),
    f_ref_(),
    f_weight_(),
    fi_ref_(),
    fi_weight_() {
}


LocalContactForceCost::~LocalContactForceCost() {
}


void LocalContactForceCost::set_f_ref(
    const std::vector<Eigen::Vector3d>& f_ref) {
  if (f_ref.size() != max_num_contacts_) {
    throw std::invalid_argument(
        "[LocalContactForceCost] invalid argument: f_ref.size() must be " 
        + std::to_string(max_num_contacts_) + "!");
  }
  f_ref_ = f_ref;
}


void LocalContactForceCost::set_f_weight(
    const std::vector<Eigen::Vector3d>& f_weight) {
  if (f_weight.size() != max_num_contacts_) {
    throw std::invalid_argument(
        "[LocalContactForceCost] invalid argument: f_weight.size() must be " 
        + std::to_string(max_num_contacts_) + "!");
  }
  f_weight_ = f_weight;
}


void LocalContactForceCost::set_fi_ref(
    const std::vector<Eigen::Vector3d>& fi_ref) {
  if (fi_ref.size() != max_num_contacts_) {
    throw std::invalid_argument(
        "[LocalContactForceCost] invalid argument: f_ref.size() must be " 
        + std::to_string(max_num_contacts_) + "!");
  }
  fi_ref_ = fi_ref;
}


void LocalContactForceCost::set_fi_weight(
    const std::vector<Eigen::Vector3d>& fi_weight) {
  if (fi_weight.size() != max_num_contacts_) {
    throw std::invalid_argument(
        "[LocalContactForceCost] invalid argument: f_weight.size() must be " 
        + std::to_string(max_num_contacts_) + "!");
  }
  fi_weight_ = fi_weight;
}


double LocalContactForceCost::evalStageCost(Robot& robot, 
                                            const ContactStatus& contact_status, 
                                            const GridInfo& grid_info,
                                            const SplitSolution& s,
                                            CostFunctionData& data) const {
  double l = 0;
  for (int i=0; i<max_num_contacts_; ++i) {
    if (contact_status.isContactActive(i)) {
      const auto& fl = s.f[i].template head<3>();
      l += (f_weight_[i].array() * (fl.array()-f_ref_[i].array()) 
                                 * (fl.array()-f_ref_[i].array())).sum();
    }
  }
  return 0.5 * l;
}


void LocalContactForceCost::evalStageCostDerivatives(
    Robot& robot, const ContactStatus& contact_status, const GridInfo& grid_info,
    const SplitSolution& s, CostFunctionData& data,
    SplitKKTResidual& kkt_residual) const {
  int dimf_stack = 0;
  for (int i=0; i<max_num_contacts_; ++i) {
    if (contact_status.isContactActive(i)) {
      const auto& fl = s.f[i].template head<3>();
      kkt_residual.lf().template segment<3>(dimf_stack).array()
          += f_weight_[i].array() * (fl.array()-f_ref_[i].array());
      switch (contact_types_[i]) {
        case ContactType::PointContact:
          dimf_stack += 3;
          break;
        case ContactType::SurfaceContact:
          dimf_stack += 6;
          break;
        default:
          break;
      }
    }
  }
}


void LocalContactForceCost::evalStageCostHessian(
    Robot& robot, const ContactStatus& contact_status, const GridInfo& grid_info,
    const SplitSolution& s, CostFunctionData& data,
    SplitKKTMatrix& kkt_matrix) const {
  int dimf_stack = 0;
  for (int i=0; i<max_num_contacts_; ++i) {
    if (contact_status.isContactActive(i)) {
      kkt_matrix.Qff().diagonal().template segment<3>(dimf_stack).noalias() 
          += f_weight_[i];
      switch (contact_types_[i]) {
        case ContactType::PointContact:
          dimf_stack += 3;
          break;
        case ContactType::SurfaceContact:
          dimf_stack += 6;
          break;
        default:
          break;
      }
    }
  }
}


double LocalContactForceCost::evalTerminalCost(Robot& robot, 
                                               const GridInfo& grid_info,  
                                               const SplitSolution& s,
                                               CostFunctionData& data) const {
  return 0;
}


void LocalContactForceCost::evalTerminalCostDerivatives(
    Robot& robot, const GridInfo& grid_info,  const SplitSolution& s,
    CostFunctionData& data, SplitKKTResidual& kkt_residual) const {
  // Do nothing.
}


void LocalContactForceCost::evalTerminalCostHessian(
    Robot& robot, const GridInfo& grid_info,  const SplitSolution& s,
    CostFunctionData& data, SplitKKTMatrix& kkt_matrix) const {
  // Do nothing.
}


double LocalContactForceCost::evalImpactCost(
    Robot& robot, const ImpactStatus& impact_status, const GridInfo& grid_info,
    const SplitSolution& s, CostFunctionData& data) const {
  double l = 0;
  for (int i=0; i<max_num_contacts_; ++i) {
    if (impact_status.isImpactActive(i)) {
      const auto& fl = s.f[i].template head<3>();
      l += (fi_weight_[i].array() * (fl.array()-fi_ref_[i].array()) 
                                  * (fl.array()-fi_ref_[i].array())).sum();
    }
  }
  return 0.5 * l;
}


void LocalContactForceCost::evalImpactCostDerivatives(
    Robot& robot, const ImpactStatus& impact_status, const GridInfo& grid_info,
    const SplitSolution& s, CostFunctionData& data, 
    SplitKKTResidual& kkt_residual) const {
  int dimf_stack = 0;
  for (int i=0; i<max_num_contacts_; ++i) {
    if (impact_status.isImpactActive(i)) {
      const auto& fl = s.f[i].template head<3>();
      kkt_residual.lf().template segment<3>(dimf_stack).array()
          += fi_weight_[i].array() * (fl.array()-fi_ref_[i].array());
      switch (contact_types_[i]) {
        case ContactType::PointContact:
          dimf_stack += 3;
          break;
        case ContactType::SurfaceContact:
          dimf_stack += 6;
          break;
        default:
          break;
      }
    }
  }
}


void LocalContactForceCost::evalImpactCostHessian(
    Robot& robot, const ImpactStatus& impact_status, const GridInfo& grid_info,
    const SplitSolution& s, CostFunctionData& data,
    SplitKKTMatrix& kkt_matrix) const {
  int dimf_stack = 0;
  for (int i=0; i<max_num_contacts_; ++i) {
    if (impact_status.isImpactActive(i)) {
      kkt_matrix.Qff().diagonal().template segment<3>(dimf_stack).noalias() 
          += fi_weight_[i];
      switch (contact_types_[i]) {
        case ContactType::PointContact:
          dimf_stack += 3;
          break;
        case ContactType::SurfaceContact:
          dimf_stack += 6;
          break;
        default:
          break;
      }
    }
  }
}

} // namespace robotoc