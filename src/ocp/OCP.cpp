#include "ocp/OCP.hpp"


namespace invdynocp {

OCP::OCP(const Robot& robot, const CostFunctionInterface* cost,
         const ConstraintsInterface* constraints, const double T, 
         const unsigned int N, const unsigned int num_proc)
  : robots_(),
    cost_(const_cast<CostFunctionInterface*>(cost)),
    constraints_(const_cast<ConstraintsInterface*>(constraints)),
    T_(T),
    N_(N),
    num_proc_(num_proc),
    dimx_(robot.dimq()+robot.dimv()),
    dimTx_(2*robot.dimv()),
    dimu_(robot.dimv()),
    dimc_(robot.dimf()+robot.dim_passive()+constraints_->dimc()) {
  for (int i=0; i<num_proc; ++i) {
    robots_.push_back(robot);
  }
}


OCP::~OCP() {
}


void OCP::solveOCP(const double t, const Eigen::VectorXd& x) {
  
}

} // namespace invdynocp