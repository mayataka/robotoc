#include "robotoc/mpc/contact_planner_base.hpp"

#include <stdexcept>
#include <iostream>
#include <cassert>


namespace robotoc {

void ContactPlannerBase::disp(std::ostream& os) const {
  Eigen::IOFormat fmt(4, 0, ", ", "\n", "", "");
  const int planning_steps = size();
  os << "Contact planner:" << std::endl;
  for (int i=0; i<planning_steps; ++i) {
    os << "  contact positions[" << i << "]: [";
    for (int j=0; j<contactPositions(i).size()-1; ++j) {
      os << "[" << contactPositions(i)[j].transpose().format(fmt) << "], ";
    }
    os << "[" << contactPositions(i)[contactPositions(i).size()-1].transpose().format(fmt) << "]";
    os << "]" << std::endl;
    os << "  contact surfaces[" << i << "]: [";
    for (int j=0; j<contactSurfaces(i).size()-1; ++j) {
      os << "[" << contactSurfaces(i)[j].row(0).format(fmt) << "]  ";
    }
    os << "[" << contactSurfaces(i)[contactSurfaces(i).size()-1].row(0).format(fmt) << "]" << std::endl;
    os << "                        ";
    for (int j=0; j<contactSurfaces(i).size()-1; ++j) {
      os << "[" << contactSurfaces(i)[j].row(1).format(fmt) << "]  ";
    }
    os << "[" << contactSurfaces(i)[contactSurfaces(i).size()-1].row(1).format(fmt) << "], " << std::endl;
    os << "                        ";
    for (int j=0; j<contactSurfaces(i).size()-1; ++j) {
      os << "[" << contactSurfaces(i)[j].row(2).format(fmt) << "], ";
    }
    os << "[" << contactSurfaces(i)[contactSurfaces(i).size()-1].row(2).format(fmt) << "]]" << std::endl;
    os << "  CoM position[" << i << "]: ["   << CoM(i).transpose().format(fmt) << "]" << std::endl;
    os << "  R[" << i << "]: [["   << R(i).row(0).format(fmt) << "]" << std::endl;
    os << "         ["             << R(i).row(1).format(fmt) << "]" << std::endl;
    os << "         ["             << R(i).row(2).format(fmt) << "]]" << std::endl;
  }
}


std::ostream& operator<<(std::ostream& os, 
                         const ContactPlannerBase& planner) {
  planner.disp(os);
  return os;
}


std::ostream& operator<<(std::ostream& os, 
                         const std::shared_ptr<ContactPlannerBase>& planner) {
  planner->disp(os);
  return os;
}

} // namespace robotoc 