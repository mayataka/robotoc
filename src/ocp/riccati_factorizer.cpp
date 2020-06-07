#include "ocp/riccati_factorizer.hpp"


namespace invdynocp {

RiccatiFactorizer::RiccatiFactorizer(const Robot& robot, 
                                     const CostFunctionInterface* cost, 
                                     const ConstraintsInterface* constraints) 
  : dim_Tx_(2*robot.dimv()),
    dim_afext_(robot.dimv()+robot.dimfmax()) {
  u_.resize(robot.dimv());
  x_res_.resize(dim_Tx_);
  lmd_res_.resize(dim_Tx_);
  afext_res_.resize(robot.dimv()+robot.dimfmax());
  C_res_.resize(robot.dim_passive()+robot.dimfmax());
  Fx_.resize(dim_Tx_, dim_Tx_);
  Fafext_.resize(dim_Tx_, dim_afext_);
  Hx_.resize(dim_Tx_, dim_Tx_);
  Hafext_.resize(dim_Tx_, dim_afext_);
  du_dq_.resize(robot.dimv(), robot.dimv());
  du_dv_.resize(robot.dimv(), robot.dimv());
  du_dafext_.resize(robot.dimv(), robot.dimv()+robot.dimfmax());
}


RiccatiFactorizer::~RiccatiFactorizer() {
}


void RiccatiFactorizer::linearizeOCP(Robot& robot, const double dtau, 
                                     const Eigen::VectorXd& x, 
                                     const Eigen::VectorXd& a, 
                                     const Eigen::VectorXd& mu, 
                                     const Eigen::VectorXd& lmd_next, 
                                     const Eigen::VectorXd& x_next) {

}


void RiccatiFactorizer::linearizeOCP(Robot& robot, const double dtau, 
                                     const Eigen::VectorXd& x, 
                                     const Eigen::VectorXd& a, 
                                     const Eigen::VectorXd& fext, 
                                     const Eigen::VectorXd& mu, 
                                     const Eigen::VectorXd& nu,
                                     const Eigen::VectorXd& lmd_next, 
                                     const Eigen::VectorXd& x_next) {
  robot.differenceState(x_next, x, x_res_);
  x_res_.head(dimv_) -= dtau * x.tail(dimv_);
  x_res_.tail(dimv_) -= dtau * a;
  if (robot.dimf() > 0) {
    robot.RNEA(x, a, fext, u_);
  }
}



} // namespace invdynocp