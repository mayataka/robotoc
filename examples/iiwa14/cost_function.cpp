#include "cost_function.hpp"


namespace idocp {
namespace iiwa14 {

CostFunction::CostFunction(const Robot* robot_ptr, const Eigen::VectorXd& q_ref)
  : CostFunctionInterface(),
    configuration_space_cost_(
        robot_ptr, q_ref, Eigen::VectorXd::Constant(robot_ptr->dimq(), 10), 
        Eigen::VectorXd::Zero(robot_ptr->dimv()),
        Eigen::VectorXd::Constant(robot_ptr->dimv(), 0.1), 
        Eigen::VectorXd::Zero(robot_ptr->dimv()),
        Eigen::VectorXd::Constant(robot_ptr->dimv(), 0.01), 
        Eigen::VectorXd::Zero(robot_ptr->dimv()),
        Eigen::VectorXd::Constant(robot_ptr->dimv(), 0.001)) {
}


CostFunction::~CostFunction() {
}


void CostFunction::lq(const Robot* robot_ptr, const double t, const double dtau,
                      const Eigen::VectorXd& q, const Eigen::VectorXd& v, 
                      const Eigen::VectorXd& a, const Eigen::VectorXd& u, 
                      Eigen::VectorXd& lq) {
  configuration_space_cost_.lq(robot_ptr, dtau, q, lq);
}


void CostFunction::lv(const Robot* robot_ptr, const double t, const double dtau,
                      const Eigen::VectorXd& q, const Eigen::VectorXd& v, 
                      const Eigen::VectorXd& a, const Eigen::VectorXd& u, 
                      Eigen::VectorXd& lv) {
  configuration_space_cost_.lv(robot_ptr, dtau, v, lv);
}


void CostFunction::la(const Robot* robot_ptr, const double t, const double dtau,
                      const Eigen::VectorXd& q, const Eigen::VectorXd& v, 
                      const Eigen::VectorXd& a, const Eigen::VectorXd& u, 
                      Eigen::VectorXd& la) {
  configuration_space_cost_.la(robot_ptr, dtau, a, la);
}


void CostFunction::lu(const Robot* robot_ptr, const double t, const double dtau,
                      const Eigen::VectorXd& q, const Eigen::VectorXd& v, 
                      const Eigen::VectorXd& a, const Eigen::VectorXd& u, 
                      Eigen::VectorXd& lu) {
  configuration_space_cost_.lu(robot_ptr, dtau, u, lu);
}


void CostFunction::lqq(const Robot* robot_ptr, const double t, 
                       const double dtau, const Eigen::VectorXd& q, 
                       const Eigen::VectorXd& v, const Eigen::VectorXd& a, 
                       const Eigen::VectorXd& u, Eigen::MatrixXd& lqq) {
  configuration_space_cost_.lqq(robot_ptr, dtau, lqq);
}


void CostFunction::lvv(const Robot* robot_ptr, const double t, 
                       const double dtau, const Eigen::VectorXd& q, 
                       const Eigen::VectorXd& v, const Eigen::VectorXd& a, 
                       const Eigen::VectorXd& u, Eigen::MatrixXd& lvv) {
  configuration_space_cost_.lvv(robot_ptr, dtau, lvv);
}

void CostFunction::laa(const Robot* robot_ptr, const double t, 
                       const double dtau, const Eigen::VectorXd& q, 
                       const Eigen::VectorXd& v, const Eigen::VectorXd& a, 
                       const Eigen::VectorXd& u, Eigen::MatrixXd& laa) {
  configuration_space_cost_.laa(robot_ptr, dtau, laa);
}


void CostFunction::luu(const Robot* robot_ptr, const double t, 
                       const double dtau, const Eigen::VectorXd& q, 
                       const Eigen::VectorXd& v, const Eigen::VectorXd& a, 
                       const Eigen::VectorXd& u, Eigen::MatrixXd& luu) {
  configuration_space_cost_.luu(robot_ptr, dtau, luu);
}


void CostFunction::phiq(const Robot* robot_ptr, const double t, 
                        const Eigen::VectorXd& q, const Eigen::VectorXd& v, 
                        Eigen::VectorXd& phiq) {
  configuration_space_cost_.phiq(robot_ptr, q, phiq);
}


void CostFunction::phiv(const Robot* robot_ptr, const double t, 
                        const Eigen::VectorXd& q, const Eigen::VectorXd& v, 
                        Eigen::VectorXd& phiv) {
  configuration_space_cost_.phiv(robot_ptr, v, phiv);
}


void CostFunction::phiqq(const Robot* robot_ptr, const double t, 
                         const Eigen::VectorXd& q, const Eigen::VectorXd& v, 
                         Eigen::MatrixXd& phiqq) {
  configuration_space_cost_.phiqq(robot_ptr, phiqq);
}


void CostFunction::phivv(const Robot* robot_ptr, const double t, 
                         const Eigen::VectorXd& q, const Eigen::VectorXd& v, 
                         Eigen::MatrixXd& phivv) {
  configuration_space_cost_.phivv(robot_ptr, phivv);
}

} // namespace iiwa14
} // namespace idocp