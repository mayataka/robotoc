#include "ocp/riccati_matrix_inverter.hpp"


namespace idocp {

RiccatiMatrixInverter::RiccatiMatrixInverter(const Robot& robot) 
  : has_floating_base_(robot.has_floating_base()),
    dimv_(robot.dimv()) {
}


void RiccatiMatrixInverter::factorize(const double dtau, 
                                    const Eigen::MatrixXd& Pqq_next, 
                                    const Eigen::MatrixXd& Pqv_next, 
                                    const Eigen::MatrixXd& Pvq_next, 
                                    const Eigen::MatrixXd& Pvv_next, 
                                    Eigen::MatrixXd& Qqq, 
                                    Eigen::MatrixXd& Qqv, 
                                    Eigen::MatrixXd& Qvq, 
                                    Eigen::MatrixXd& Qvv) {
  if (has_floating_base_) {

  }
  else {
    Qqq.noalias() += Pqq_next;
    Qqv.noalias() += dtau * Pqq_next;
    Qqv.noalias() += Pqv_next;
    Qvq.noalias() += dtau * Pqq_next;
    Qvq.noalias() += Pvq_next;
    Qvv.noalias() += (dtau*dtau) * Pqq_next;
    Qvv.noalias() += dtau * (Pqv_next + Pvq_next);
    Qvv.noalias() += Pvv_next;
  }
}


void RiccatiMatrixFactorizer::factorize(const double dtau, 
                                        const Eigen::MatrixXd& Pqv_next, 
                                        const Eigen::MatrixXd& Pvv_next, 
                                        Eigen::MatrixXd& Qqa, 
                                        Eigen::MatrixXd& Qva) {
  if (has_floating_base_) {

  }
  else {
    Qqa.noalias() += dtau * Pqv_next;
    Qva.noalias() += (dtau*dtau) * Pqv_next;
    Qva.noalias() += dtau * Pvv_next;
  }
}


void RiccatiMatrixFactorizer::factorize(const double dtau, 
                                        const Eigen::MatrixXd& Pvv_next, 
                                        Eigen::MatrixXd& Qaa) {
  Qaa_.noalias() += (dtau*dtau) * Pvv_next;
}

} // namespace idocp