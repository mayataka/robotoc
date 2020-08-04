#include "idocp/ocp/riccati_matrix_factorizer.hpp"

#include <assert.h>


namespace idocp {

RiccatiMatrixFactorizer::RiccatiMatrixFactorizer(const Robot& robot) 
  : has_floating_base_(robot.has_floating_base()),
    dimv_(robot.dimv()),
    dintegrate_dq_(),
    dintegrate_dv_() {
  if (robot.has_floating_base()) {
    dintegrate_dq_.resize(robot.dimv(), robot.dimv());
    dintegrate_dv_.resize(robot.dimv(), robot.dimv());
  }
}


RiccatiMatrixFactorizer::RiccatiMatrixFactorizer() 
  : has_floating_base_(false),
    dimv_(0),
    dintegrate_dq_(),
    dintegrate_dv_() {
}


RiccatiMatrixFactorizer::~RiccatiMatrixFactorizer() {
}


void RiccatiMatrixFactorizer::setIntegrationSensitivities(
    const Robot& robot, const double dtau, const Eigen::VectorXd& q,
    const Eigen::VectorXd& v) {
  assert(dtau > 0);
  assert(q.size() == robot.dimq());
  assert(v.size() == robot.dimv());
  if (has_floating_base_) {
    robot.dIntegrateConfiguration(q, v, dtau, dintegrate_dq_, dintegrate_dv_);
  }
}


void RiccatiMatrixFactorizer::factorize(const double dtau, 
                                        const Eigen::MatrixXd& Pqq_next, 
                                        const Eigen::MatrixXd& Pqv_next, 
                                        const Eigen::MatrixXd& Pvq_next, 
                                        const Eigen::MatrixXd& Pvv_next, 
                                        Eigen::MatrixXd& Qqq, 
                                        Eigen::MatrixXd& Qqv, 
                                        Eigen::MatrixXd& Qvq, 
                                        Eigen::MatrixXd& Qvv) {
  assert(dtau > 0);
  assert(Pqq_next.rows() == dimv_);
  assert(Pqq_next.cols() == dimv_);
  assert(Pqv_next.rows() == dimv_);
  assert(Pqv_next.cols() == dimv_);
  assert(Pvq_next.rows() == dimv_);
  assert(Pvq_next.cols() == dimv_);
  assert(Pvv_next.rows() == dimv_);
  assert(Pvv_next.cols() == dimv_);
  assert(Qqq.rows() == dimv_);
  assert(Qqq.cols() == dimv_);
  assert(Qqv.rows() == dimv_);
  assert(Qqv.cols() == dimv_);
  assert(Qvq.rows() == dimv_);
  assert(Qvq.cols() == dimv_);
  assert(Qvv.rows() == dimv_);
  assert(Qvv.cols() == dimv_);
  if (has_floating_base_) {
    Qqq.topLeftCorner<6, 6>().noalias() 
        += dintegrate_dq_.topLeftCorner<6, 6>().transpose() 
            * Pqq_next.topLeftCorner<6, 6>()
            * dintegrate_dq_.topLeftCorner<6, 6>();
    Qqq.topRightCorner(6, dimv_-6).noalias()
        += dintegrate_dq_.topLeftCorner<6, 6>().transpose()
            * Pqq_next.topRightCorner(6, dimv_-6);
    Qqq.bottomLeftCorner(dimv_-6, 6).noalias()
        += Pqq_next.bottomLeftCorner(dimv_-6, 6) 
            * dintegrate_dq_.topLeftCorner<6, 6>();
    Qqq.bottomRightCorner(dimv_-6, dimv_-6).noalias()
        += Pqq_next.bottomRightCorner(dimv_-6, dimv_-6);

    Qqv.topLeftCorner<6, 6>().noalias() 
        += dtau * dintegrate_dq_.topLeftCorner<6, 6>().transpose() 
                * Pqq_next.topLeftCorner<6, 6>() 
                * dintegrate_dv_.topLeftCorner<6, 6>();
    Qqv.topRightCorner(6, dimv_-6).noalias()
        += dtau * dintegrate_dq_.topLeftCorner<6, 6>().transpose()
                * Pqq_next.topRightCorner(6, dimv_-6);
    Qqv.bottomLeftCorner(dimv_-6, 6).noalias()
        += dtau * Pqq_next.bottomLeftCorner(dimv_-6, 6) 
                * dintegrate_dv_.topLeftCorner<6, 6>();
    Qqv.bottomRightCorner(dimv_-6, dimv_-6).noalias()
        += dtau * Pqq_next.bottomRightCorner(dimv_-6, dimv_-6);
    Qqv.topRows<6>().noalias() 
        += dintegrate_dq_.topLeftCorner<6, 6>().transpose()
            * Pqv_next.topRows<6>();
    Qqv.bottomRows(dimv_-6).noalias() += Pqv_next.bottomRows(dimv_-6);

    Qvq.topLeftCorner<6, 6>().noalias() 
        += dtau * dintegrate_dv_.topLeftCorner<6, 6>().transpose() 
                * Pqq_next.topLeftCorner<6, 6>() 
                * dintegrate_dq_.topLeftCorner<6, 6>();
    Qvq.topRightCorner(6, dimv_-6).noalias()
        += dtau *dintegrate_dv_.topLeftCorner<6, 6>().transpose()
                * Pqq_next.topRightCorner(6, dimv_-6);
    Qvq.bottomLeftCorner(dimv_-6, 6).noalias()
        += dtau * Pqq_next.bottomLeftCorner(dimv_-6, 6) 
                * dintegrate_dq_.topLeftCorner<6, 6>();
    Qvq.bottomRightCorner(dimv_-6, dimv_-6).noalias()
        += dtau * Pqq_next.bottomRightCorner(dimv_-6, dimv_-6);
    Qvq.leftCols<6>().noalias() 
        += Pvq_next.leftCols<6>() * dintegrate_dq_.topLeftCorner<6, 6>();
    Qvq.rightCols(dimv_-6).noalias() += Pvq_next.rightCols(dimv_-6);

    Qvv.topLeftCorner<6, 6>().noalias() 
        += (dtau*dtau) * dintegrate_dv_.topLeftCorner<6, 6>().transpose() 
                       * Pqq_next.topLeftCorner<6, 6>()
                       * dintegrate_dv_.topLeftCorner<6, 6>();
    Qvv.topRightCorner(6, dimv_-6).noalias()
        += (dtau*dtau) * dintegrate_dv_.topLeftCorner<6, 6>().transpose()
                       * Pqq_next.topRightCorner(6, dimv_-6);
    Qvv.bottomLeftCorner(dimv_-6, 6).noalias()
        += (dtau*dtau) * Pqq_next.bottomLeftCorner(dimv_-6, 6) 
                       * dintegrate_dv_.topLeftCorner<6, 6>();
    Qvv.bottomRightCorner(dimv_-6, dimv_-6).noalias()
        += (dtau*dtau) * Pqq_next.bottomRightCorner(dimv_-6, dimv_-6);
    Qvv.leftCols<6>().noalias() 
        += dtau * Pvq_next.leftCols<6>() * dintegrate_dv_.topLeftCorner<6, 6>();
    Qvv.rightCols(dimv_-6).noalias() += dtau * Pvq_next.rightCols(dimv_-6);
    Qvv.topRows<6>().noalias() 
        += dtau * dintegrate_dv_.topLeftCorner<6, 6>().transpose()
                * Pqv_next.topRows<6>();
    Qvv.bottomRows(dimv_-6).noalias() += dtau * Pqv_next.bottomRows(dimv_-6);
    Qvv.noalias() += Pvv_next;
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
  assert(dtau > 0);
  assert(Pqv_next.rows() == dimv_);
  assert(Pqv_next.cols() == dimv_);
  assert(Pvv_next.rows() == dimv_);
  assert(Pvv_next.cols() == dimv_);
  assert(Qqa.rows() == dimv_);
  assert(Qqa.cols() == dimv_);
  assert(Qva.rows() == dimv_);
  assert(Qva.cols() == dimv_);
  if (has_floating_base_) {
    Qqa.topRows<6>().noalias() 
        += dtau * dintegrate_dq_.topLeftCorner<6, 6>().transpose()
                * Pqv_next.topRows<6>();
    Qqa.bottomRows(dimv_-6).noalias() += dtau * Pqv_next.bottomRows(dimv_-6);
    Qva.topRows<6>().noalias() 
        += (dtau*dtau) * dintegrate_dv_.topLeftCorner<6, 6>().transpose()
                       * Pqv_next.topRows<6>();
    Qva.bottomRows(dimv_-6).noalias() 
        += (dtau*dtau) * Pqv_next.bottomRows(dimv_-6);
    Qva.noalias() += dtau * Pvv_next;
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
  assert(dtau > 0);
  assert(Pvv_next.rows() == dimv_);
  assert(Pvv_next.cols() == dimv_);
  assert(Qaa.rows() == dimv_);
  assert(Qaa.cols() == dimv_);
  Qaa.noalias() += (dtau*dtau) * Pvv_next;
}

} // namespace idocp