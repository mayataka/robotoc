#include "ocp/riccati_matrix_inverter.hpp"

#include <assert.h>
#include "Eigen/LU"


namespace idocp {

RiccatiMatrixInverter::RiccatiMatrixInverter(const Robot& robot) 
  : has_floating_base_(robot.has_floating_base()),
    dimv_(robot.dimv()),
    max_dimf_(robot.max_dimf()),
    dimf_(0),
    dim_passive_(robot.dim_passive()),
    Sff_inv_(Eigen::MatrixXd::Zero(robot.max_dimf(), robot.max_dimf())),
    Saa_(Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv())),
    Saa_inv_(Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv())),
    Saf_(Eigen::MatrixXd::Zero(robot.dimv(), robot.max_dimf())),
    Sac_(Eigen::MatrixXd::Zero(robot.dimv(), 
                               robot.max_dimf()+robot.dim_passive())),
    Sfc_(Eigen::MatrixXd::Zero(robot.max_dimf(), 
                               robot.max_dimf()+robot.dim_passive())),
    Scc_(Eigen::MatrixXd::Zero(robot.max_dimf()+robot.dim_passive(), 
                               robot.max_dimf()+robot.dim_passive())),
    Scc_inv_(Eigen::MatrixXd::Zero(robot.max_dimf()+robot.dim_passive(), 
                                   robot.max_dimf()+robot.dim_passive())) {
  if (robot.max_dimf() == 0) {
    Saa_.resize(0, 0);
  }
}


RiccatiMatrixInverter::~RiccatiMatrixInverter() {
}


void RiccatiMatrixInverter::setContactStatus(const Robot& robot) {
  dimf_ = robot.dimf();
}


void RiccatiMatrixInverter::precompute(const Eigen::MatrixXd& Qff, 
                                       const Eigen::MatrixXd& Qaf) {
  assert(Qff.rows() == max_dimf_);
  assert(Qff.cols() == max_dimf_);
  assert(Qaf.rows() == dimv_);
  assert(Qaf.cols() == max_dimf_);
  if (dimf_ > 0) {
    Sff_inv_.topLeftCorner(dimf_, dimf_) 
        = Qff.topLeftCorner(dimf_, dimf_)
             .llt().solve(Eigen::MatrixXd::Identity(dimf_, dimf_));
    Saa_ = - Qaf.leftCols(dimf_) * Sff_inv_.topLeftCorner(dimf_, dimf_) 
                                 * Qaf.leftCols(dimf_).transpose();
    Sac_.leftCols(dimf_) 
        = Qaf.leftCols(dimf_) * Sff_inv_.topLeftCorner(dimf_, dimf_);
  }
}


void RiccatiMatrixInverter::invert(const Eigen::MatrixXd& Qqa, 
                                   const Eigen::MatrixXd& Qva,
                                   const Eigen::MatrixXd& Qaa, 
                                   const Eigen::VectorXd& la,
                                   Eigen::MatrixXd& Kaq, Eigen::MatrixXd& Kav, 
                                   Eigen::VectorXd& ka) {
  assert(!has_floating_base_);
  assert(dim_passive_ == 0);
  assert(Qqa.rows() == dimv_);
  assert(Qqa.cols() == dimv_);
  assert(Qva.rows() == dimv_);
  assert(Qva.cols() == dimv_);
  assert(Qaa.rows() == dimv_);
  assert(Qaa.cols() == dimv_);
  assert(la.size() == dimv_);
  assert(Kaq.rows() == dimv_);
  assert(Kaq.cols() == dimv_);
  assert(Kav.rows() == dimv_);
  assert(Kav.cols() == dimv_);
  assert(ka.size() == dimv_);
  assert(dimf_ == 0);
  Saa_inv_ = Qaa.llt().solve(Eigen::MatrixXd::Identity(dimv_, dimv_));
  Kaq = - Saa_inv_ * Qqa.transpose();
  Kav = - Saa_inv_ * Qva.transpose();
  ka = - Saa_inv_ * la;
}


void RiccatiMatrixInverter::invert(const Eigen::MatrixXd& Qqa, 
                                   const Eigen::MatrixXd& Qva,
                                   const Eigen::MatrixXd& Qaa, 
                                   const Eigen::MatrixXd& Qqf,
                                   const Eigen::MatrixXd& Qvf, 
                                   const Eigen::MatrixXd& Cq,
                                   const Eigen::MatrixXd& Cv,
                                   const Eigen::MatrixXd& Ca,
                                   const Eigen::VectorXd& la,
                                   const Eigen::VectorXd& lf, 
                                   const Eigen::VectorXd& C_res,
                                   Eigen::MatrixXd& Kaq, Eigen::MatrixXd& Kav, 
                                   Eigen::MatrixXd& Kfq, Eigen::MatrixXd& Kfv, 
                                   Eigen::MatrixXd& Kmuq, Eigen::MatrixXd& Kmuv, 
                                   Eigen::VectorXd& ka, Eigen::VectorXd& kf,
                                   Eigen::VectorXd& kmu) {
  assert(!has_floating_base_);
  assert(dim_passive_ == 0);
  assert(Qqa.rows() == dimv_);
  assert(Qqa.cols() == dimv_);
  assert(Qva.rows() == dimv_);
  assert(Qva.cols() == dimv_);
  assert(Qaa.rows() == dimv_);
  assert(Qaa.cols() == dimv_);
  assert(Qqf.rows() == dimv_);
  assert(Qqf.cols() == max_dimf_);
  assert(Qvf.rows() == dimv_);
  assert(Qvf.cols() == max_dimf_);
  assert(la.size() == dimv_);
  assert(lf.size() == max_dimf_);
  assert(C_res.size() == max_dimf_+dim_passive_);
  assert(Kaq.rows() == dimv_);
  assert(Kaq.cols() == dimv_);
  assert(Kav.rows() == dimv_);
  assert(Kav.cols() == dimv_);
  assert(Kfq.rows() == max_dimf_);
  assert(Kfq.cols() == dimv_);
  assert(Kfv.rows() == max_dimf_);
  assert(Kfv.cols() == dimv_);
  assert(Kmuq.rows() == max_dimf_+dim_passive_);
  assert(Kmuq.cols() == dimv_);
  assert(Kmuv.rows() == max_dimf_+dim_passive_);
  assert(Kmuv.cols() == dimv_);
  assert(ka.size() == dimv_);
  assert(kf.size() == max_dimf_);
  assert(kmu.size() == max_dimf_+dim_passive_);
  // Inverts the coefficient matrix of the decision variables.
  Saa_.noalias() += Qaa;
  Saa_inv_ = Saa_.llt().solve(Eigen::MatrixXd::Identity(dimv_, dimv_));
  Saf_.leftCols(dimf_) = Saa_inv_ * Sac_.leftCols(dimf_);
  Scc_.topLeftCorner(dimf_, dimf_) 
      = - Ca.topRows(dimf_) * Saa_inv_ * Ca.topRows(dimf_).transpose();
  Scc_inv_.topLeftCorner(dimf_, dimf_) 
      = Scc_.topLeftCorner(dimf_, dimf_)
            .llt().solve(Eigen::MatrixXd::Identity(dimf_, dimf_));
  Sac_.leftCols(dimf_) 
      = Saa_inv_ * Ca.topRows(dimf_).transpose() 
                 * Scc_inv_.topLeftCorner(dimf_, dimf_);
  Sfc_.topLeftCorner(dimf_, dimf_) 
      = - Saf_.leftCols(dimf_).transpose() 
          * Ca.topRows(dimf_).transpose() 
          * Scc_inv_.topLeftCorner(dimf_, dimf_);
  Saa_inv_.noalias() 
      += Sac_.leftCols(dimf_) * Scc_.topLeftCorner(dimf_, dimf_)
                              * Sac_.leftCols(dimf_).transpose();
  Sff_inv_.topLeftCorner(dimf_, dimf_).noalias() 
      += Saf_.leftCols(dimf_).transpose() * Saa_ * Saf_.leftCols(dimf_);
  Sff_inv_.topLeftCorner(dimf_, dimf_).noalias() 
      += Sfc_.topLeftCorner(dimf_, dimf_) 
          * Scc_.topLeftCorner(dimf_, dimf_) 
          * Sfc_.topLeftCorner(dimf_, dimf_).transpose();
  Saf_.leftCols(dimf_).noalias() 
      -= Sac_.leftCols(dimf_) * Scc_.topLeftCorner(dimf_, dimf_)
                              * Sfc_.topLeftCorner(dimf_, dimf_).transpose();
  // Computes the state feedback gain.
  Kaq = - Saa_inv_ * Qqa.transpose();
  Kaq.noalias() += Saf_.leftCols(dimf_) * Qqf.leftCols(dimf_).transpose();
  Kaq.noalias() += Sac_.leftCols(dimf_) * Cq.topRows(dimf_);
  Kav = - Saa_inv_ * Qva.transpose();
  Kav.noalias() += Saf_.leftCols(dimf_) * Qvf.leftCols(dimf_).transpose();
  Kav.noalias() += Sac_.leftCols(dimf_) * Cv.topRows(dimf_);
  Kfq.topRows(dimf_) = Saf_.leftCols(dimf_).transpose() * Qqa.transpose();
  Kfq.topRows(dimf_).noalias() 
      -= Sff_inv_.topLeftCorner(dimf_, dimf_) 
          * Qqf.leftCols(dimf_).transpose();
  Kfq.topRows(dimf_).noalias() 
      += Sfc_.topLeftCorner(dimf_, dimf_) * Cq.topRows(dimf_);
  Kfv.topRows(dimf_) = Saf_.leftCols(dimf_).transpose() * Qva.transpose();
  Kfv.topRows(dimf_).noalias() 
      -= Sff_inv_.topLeftCorner(dimf_, dimf_) 
          * Qvf.leftCols(dimf_).transpose();
  Kfv.topRows(dimf_).noalias() 
      += Sfc_.topLeftCorner(dimf_, dimf_) * Cv.topRows(dimf_);
  Kmuq.topRows(dimf_) = Sac_.leftCols(dimf_).transpose() * Qqa.transpose();
  Kmuq.topRows(dimf_).noalias() 
      += Sfc_.topLeftCorner(dimf_, dimf_).transpose() 
          * Qqf.leftCols(dimf_).transpose();
  Kmuq.topRows(dimf_).noalias() 
      -= Scc_.topLeftCorner(dimf_, dimf_) * Cq.topRows(dimf_);
  Kmuv.topRows(dimf_) = Sac_.leftCols(dimf_).transpose() * Qva.transpose();
  Kmuv.topRows(dimf_).noalias() 
      += Sfc_.topLeftCorner(dimf_, dimf_).transpose() 
          * Qvf.leftCols(dimf_).transpose();
  Kmuv.topRows(dimf_).noalias() 
      -= Scc_.topLeftCorner(dimf_, dimf_) * Cv.topRows(dimf_);
  // Computes the state feedforward terms.
  ka = - Saa_inv_ * la;
  ka.noalias() += Saf_.leftCols(dimf_) * lf.head(dimf_);
  ka.noalias() += Sac_.leftCols(dimf_) * C_res.head(dimf_);
  kf.topRows(dimf_) = Saf_.leftCols(dimf_).transpose() * la;
  kf.topRows(dimf_).noalias() 
      -= Sff_inv_.topLeftCorner(dimf_, dimf_) * lf.head(dimf_);
  kf.topRows(dimf_).noalias() 
      += Sfc_.topLeftCorner(dimf_, dimf_) * C_res.head(dimf_);
  kmu.topRows(dimf_) = Sac_.leftCols(dimf_).transpose() * la;
  kmu.topRows(dimf_).noalias() 
      += Sfc_.topLeftCorner(dimf_, dimf_).transpose() * lf.head(dimf_);
  kmu.topRows(dimf_).noalias() 
      -= Scc_.topLeftCorner(dimf_, dimf_) * C_res.head(dimf_);
}

} // namespace idocp