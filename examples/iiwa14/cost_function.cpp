#include "cost_function.hpp"


namespace idocp {
namespace iiwa14 {

CostFunction::CostFunction(const Robot& robot, const Eigen::VectorXd& q_ref)
  : CostFunctionInterface(),
    joint_space_cost_(robot, Eigen::VectorXd::Constant(robot.dimq(), 10), 
                      Eigen::VectorXd::Constant(robot.dimv(), 1), 
                      Eigen::VectorXd::Constant(robot.dimv(), 0.01), 
                      Eigen::VectorXd::Constant(robot.dimv(), 0.0), 
                      Eigen::VectorXd::Constant(robot.dimq(), 10), 
                      Eigen::VectorXd::Constant(robot.dimv(), 1)) {
  joint_space_cost_.set_qref(q_ref);
}


CostFunction::~CostFunction() {
}

double CostFunction::l(const Robot& robot, const double t, const double dtau,
                       const Eigen::VectorXd& q, const Eigen::VectorXd& v, 
                       const Eigen::VectorXd& a, const Eigen::VectorXd& u) {
  double l = 0;
  l += joint_space_cost_.l(robot, dtau, q, v, a, u);
  return l;
}


double CostFunction::phi(const Robot& robot, const double t, 
                         const Eigen::VectorXd& q, const Eigen::VectorXd& v) {
  double phi = 0;
  phi += joint_space_cost_.phi(robot, q, v);
  return phi;
}


void CostFunction::lq(const Robot& robot, const double t, const double dtau,
                      const Eigen::VectorXd& q, const Eigen::VectorXd& v, 
                      const Eigen::VectorXd& a, Eigen::VectorXd& lq) {
  joint_space_cost_.lq(robot, dtau, q, lq);
}


void CostFunction::lv(const Robot& robot, const double t, const double dtau,
                      const Eigen::VectorXd& q, const Eigen::VectorXd& v, 
                      const Eigen::VectorXd& a, Eigen::VectorXd& lv) {
  joint_space_cost_.lv(robot, dtau, v, lv);
}


void CostFunction::la(const Robot& robot, const double t, const double dtau,
                      const Eigen::VectorXd& q, const Eigen::VectorXd& v, 
                      const Eigen::VectorXd& a, Eigen::VectorXd& la) {
  joint_space_cost_.la(robot, dtau, a, la);
}


void CostFunction::lu(const Robot& robot, const double t, const double dtau,
                      const Eigen::VectorXd& u, Eigen::VectorXd& lu) {
  joint_space_cost_.lu(robot, dtau, u, lu);
}


void CostFunction::lqq(const Robot& robot, const double t, const double dtau, 
                       const Eigen::VectorXd& q, const Eigen::VectorXd& v, 
                       const Eigen::VectorXd& a, Eigen::MatrixXd& lqq) {
  joint_space_cost_.lqq(robot, dtau, lqq);
}


void CostFunction::lvv(const Robot& robot, const double t, const double dtau, 
                       const Eigen::VectorXd& q, const Eigen::VectorXd& v, 
                       const Eigen::VectorXd& a, Eigen::MatrixXd& lvv) {
  joint_space_cost_.lvv(robot, dtau, lvv);
}


void CostFunction::laa(const Robot& robot, const double t, const double dtau, 
                       const Eigen::VectorXd& q, const Eigen::VectorXd& v, 
                       const Eigen::VectorXd& a, Eigen::MatrixXd& laa) {
  joint_space_cost_.laa(robot, dtau, laa);
}


void CostFunction::luu(const Robot& robot, const double t, const double dtau, 
                       const Eigen::VectorXd& u, Eigen::MatrixXd& luu) {
  joint_space_cost_.luu(robot, dtau, luu);
}


void CostFunction::phiq(const Robot& robot, const double t, 
                        const Eigen::VectorXd& q, const Eigen::VectorXd& v, 
                        Eigen::VectorXd& phiq) {
  joint_space_cost_.phiq(robot, q, phiq);
}


void CostFunction::phiv(const Robot& robot, const double t, 
                        const Eigen::VectorXd& q, const Eigen::VectorXd& v, 
                        Eigen::VectorXd& phiv) {
  joint_space_cost_.phiv(robot, v, phiv);
}


void CostFunction::phiqq(const Robot& robot, const double t, 
                         const Eigen::VectorXd& q, const Eigen::VectorXd& v, 
                         Eigen::MatrixXd& phiqq) {
  joint_space_cost_.phiqq(robot, phiqq);
}


void CostFunction::phivv(const Robot& robot, const double t, 
                         const Eigen::VectorXd& q, const Eigen::VectorXd& v, 
                         Eigen::MatrixXd& phivv) {
  joint_space_cost_.phivv(robot, phivv);
}

} // namespace iiwa14
} // namespace idocp