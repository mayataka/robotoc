#include "cost/configuration_space_cost.hpp"

#include <assert.h>


namespace invdynocp {

ConfigurationSpaceCost::ConfigurationSpaceCost(const Robot* robot_ptr, 
                                               const Eigen::VectorXd& q_weight,  
                                               const Eigen::VectorXd& v_weight,  
                                               const Eigen::VectorXd& a_weight,  
                                               const Eigen::VectorXd& u_weight) 
  : dimq_(robot_ptr->dimq()),
    dimv_(robot_ptr->dimv()),
    q_ref_(Eigen::VectorXd::Zero(robot_ptr->dimq())),
    v_ref_(Eigen::VectorXd::Zero(robot_ptr->dimv())),
    a_ref_(Eigen::VectorXd::Zero(robot_ptr->dimv())),
    u_ref_(Eigen::VectorXd::Zero(robot_ptr->dimv())),
    q_weight_(q_weight),
    v_weight_(v_weight),
    a_weight_(a_weight),
    u_weight_(u_weight) {
}


ConfigurationSpaceCost::ConfigurationSpaceCost(const Robot* robot_ptr, 
                                               const Eigen::VectorXd& q_ref, 
                                               const Eigen::VectorXd& q_weight, 
                                               const Eigen::VectorXd& v_ref, 
                                               const Eigen::VectorXd& v_weight, 
                                               const Eigen::VectorXd& a_ref, 
                                               const Eigen::VectorXd& a_weight, 
                                               const Eigen::VectorXd& u_ref, 
                                               const Eigen::VectorXd& u_weight) 
  : dimq_(robot_ptr->dimq()),
    dimv_(robot_ptr->dimv()),
    q_ref_(q_ref),
    q_weight_(q_weight),
    v_ref_(v_ref),
    v_weight_(v_weight),
    a_ref_(a_ref),
    a_weight_(a_weight),
    u_ref_(u_ref),
    u_weight_(u_weight) {
}


void ConfigurationSpaceCost::lq(const Robot* robot_ptr, const double dtau, 
                                const Eigen::VectorXd& q, Eigen::VectorXd& lq) {
  for (int i=0; i<dimq_; ++i) {
    lq.coeffRef(i) = dtau * q_weight_.coeff(i) * (q.coeff(i)-q_ref_.coeff(i));
  }
}


void ConfigurationSpaceCost::lv(const Robot* robot_ptr, const double dtau, 
                                const Eigen::VectorXd& v, Eigen::VectorXd& lv) {
  for (int i=0; i<dimv_; ++i) {
    lv.coeffRef(i) = dtau * v_weight_.coeff(i) * (v.coeff(i)-v_ref_.coeff(i));
  }
}


void ConfigurationSpaceCost::la(const Robot* robot_ptr, const double dtau, 
                                const Eigen::VectorXd& a, Eigen::VectorXd& la) {
  for (int i=0; i<dimv_; ++i) {
    la.coeffRef(i) = dtau * a_weight_.coeff(i) * (a.coeff(i)-a_ref_.coeff(i));
  }
}


void ConfigurationSpaceCost::lu(const Robot* robot_ptr, const double dtau, 
                                const Eigen::VectorXd& u, Eigen::VectorXd& lu) {
  for (int i=0; i<dimv_; ++i) {
    lu.coeffRef(i) = dtau * u_weight_.coeff(i) * (u.coeff(i) - u_ref_.coeff(i));
  }
}


void ConfigurationSpaceCost::lqq(const Robot* robot_ptr, const double dtau, 
                                 Eigen::MatrixXd& lqq) {
  for (int i=0; i<dimq_; ++i) {
    lqq.coeffRef(i, i) += dtau * q_weight_.coeff(i);
  }
}


void ConfigurationSpaceCost::lvv(const Robot* robot_ptr, const double dtau, 
                                 Eigen::MatrixXd& lvv) {
  for (int i=0; i<dimv_; ++i) {
    lvv.coeffRef(i, i) += dtau * v_weight_.coeff(i);
  }
}


void ConfigurationSpaceCost::laa(const Robot* robot_ptr, const double dtau, 
                                 Eigen::MatrixXd& laa) {
  for (int i=0; i<dimv_; ++i) {
    laa.coeffRef(i, i) += dtau * a_weight_.coeff(i);
  }
}

void ConfigurationSpaceCost::luu(const Robot* robot_ptr, const double dtau, 
                                 Eigen::MatrixXd& luu) {
  for (int i=0; i<dimv_; ++i) {
    luu.coeffRef(i, i) = dtau * u_weight_.coeff(i);
  }
}


void ConfigurationSpaceCost::phiq(const Robot* robot_ptr, 
                                  const Eigen::VectorXd& q, 
                                  Eigen::VectorXd& phiq) {
  for (int i=0; i<dimq_; ++i) {
    phiq.coeffRef(i) = q_weight_.coeff(i) * (q.coeff(i)-q_ref_.coeff(i));
  }
}


void ConfigurationSpaceCost::phiv(const Robot* robot_ptr, 
                                  const Eigen::VectorXd& v, 
                                  Eigen::VectorXd& phiv) {
  for (int i=0; i<dimv_; ++i) {
    phiv.coeffRef(i) = v_weight_.coeff(i) * (v.coeff(i)-v_ref_.coeff(i));
  }
}


void ConfigurationSpaceCost::phiqq(const Robot* robot_ptr, 
                                   Eigen::MatrixXd& phiqq) {
  for (int i=0; i<dimq_; ++i) {
    phiqq.coeffRef(i, i) = q_weight_.coeff(i);
  }
}


void ConfigurationSpaceCost::phivv(const Robot* robot_ptr, 
                                   Eigen::MatrixXd& phivv) {
  for (int i=0; i<dimv_; ++i) {
    phivv.coeffRef(i, i) = v_weight_.coeff(i);
  }
}


} // namespace invdynocp