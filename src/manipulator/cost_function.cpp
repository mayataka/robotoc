#include "manipulator/cost_function.hpp"

#include <assert.h>


namespace idocp {
namespace manipulator {

CostFunction::CostFunction(const Robot& robot)
  : CostFunctionInterface(),
    joint_space_cost_(robot, Eigen::VectorXd::Constant(robot.dimq(), 10), 
                      Eigen::VectorXd::Constant(robot.dimv(), 1), 
                      Eigen::VectorXd::Constant(robot.dimv(), 0.01), 
                      Eigen::VectorXd::Constant(robot.dimv(), 0.0), 
                      Eigen::VectorXd::Constant(robot.dimq(), 10), 
                      Eigen::VectorXd::Constant(robot.dimv(), 1)) {
}


CostFunction::~CostFunction() {
}


void CostFunction::set_q_ref(const Eigen::VectorXd& q_ref) {
  joint_space_cost_.set_q_ref(q_ref);
}


void CostFunction::set_v_ref(const Eigen::VectorXd& v_ref) {
  joint_space_cost_.set_v_ref(v_ref);
}


void CostFunction::set_a_ref(const Eigen::VectorXd& a_ref) {
  joint_space_cost_.set_a_ref(a_ref);
}


void CostFunction::set_u_ref(const Eigen::VectorXd& q_ref) {
  joint_space_cost_.set_u_ref(q_ref);
}


void CostFunction::set_q_weight(const Eigen::VectorXd& q_weight) {
  joint_space_cost_.set_q_weight(q_weight);
}


void CostFunction::set_v_weight(const Eigen::VectorXd& v_weight) {
  joint_space_cost_.set_v_weight(v_weight);
}


void CostFunction::set_a_weight(const Eigen::VectorXd& a_weight) {
  joint_space_cost_.set_a_weight(a_weight);
}


void CostFunction::set_u_weight(const Eigen::VectorXd& u_weight) {
  joint_space_cost_.set_u_weight(u_weight);
}


void CostFunction::set_qf_weight(const Eigen::VectorXd& qf_weight) {
  joint_space_cost_.set_qf_weight(qf_weight);
}


void CostFunction::set_vf_weight(const Eigen::VectorXd& vf_weight) {
  joint_space_cost_.set_vf_weight(vf_weight);
}


double CostFunction::l(const Robot& robot, const double t, 
                                  const double dtau, const Eigen::VectorXd& q, 
                                  const Eigen::VectorXd& v, 
                                  const Eigen::VectorXd& a, 
                                  const Eigen::VectorXd& u) {
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


} // namespace manipulator
} // namespace idocp