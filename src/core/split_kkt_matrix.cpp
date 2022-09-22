#include "robotoc/core/split_kkt_matrix.hpp"

#include <random>

namespace robotoc {

SplitKKTMatrix::SplitKKTMatrix(const Robot& robot) 
  : Fxx(Eigen::MatrixXd::Zero(2*robot.dimv(), 2*robot.dimv())),
    Fvu(Eigen::MatrixXd::Zero(robot.dimv(), robot.dimu())),
    Qxx(Eigen::MatrixXd::Zero(2*robot.dimv(), 2*robot.dimv())),
    Qaa(Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv())),
    Qdvdv(Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv())),
    Qxu(Eigen::MatrixXd::Zero(2*robot.dimv(), robot.dimu())),
    Quu(Eigen::MatrixXd::Zero(robot.dimu(), robot.dimu())),
    fx(Eigen::VectorXd::Zero(2*robot.dimv())),
    Qtt(0),
    Qtt_prev(0),
    hx(Eigen::VectorXd::Zero(2*robot.dimv())),
    ha(Eigen::VectorXd::Zero(robot.dimv())),
    hu(Eigen::VectorXd::Zero(robot.dimu())),
    Qff_full_(Eigen::MatrixXd::Zero(robot.max_dimf(), robot.max_dimf())),
    Qqf_full_(Eigen::MatrixXd::Zero(robot.dimv(), robot.max_dimf())),
    hf_full_(Eigen::VectorXd::Zero(robot.max_dimf())),
    has_floating_base_(robot.hasFloatingBase()),
    dimv_(robot.dimv()), 
    dimx_(2*robot.dimv()), 
    dimu_(robot.dimu()), 
    dimf_(0) {
}


SplitKKTMatrix::SplitKKTMatrix() 
  : Fxx(),
    Fvu(),
    Qxx(),
    Qaa(),
    Qdvdv(),
    Qxu(),
    Quu(),
    fx(),
    Qtt(0),
    Qtt_prev(0),
    hx(),
    ha(),
    hu(),
    Qff_full_(),
    Qqf_full_(),
    hf_full_(),
    has_floating_base_(false),
    dimv_(0), 
    dimx_(0), 
    dimu_(0), 
    dimf_(0) {
}


bool SplitKKTMatrix::isDimensionConsistent() const {
  if (Fxx.rows() != 2*dimv_) return false;
  if (Fxx.cols() != 2*dimv_) return false;
  if (Fvu.rows() != dimv_) return false;
  if (Fvu.cols() != dimu_) return false;
  if (Qxx.rows() != 2*dimv_) return false;
  if (Qxx.cols() != 2*dimv_) return false;
  if (Qaa.rows() != dimv_) return false;
  if (Qaa.cols() != dimv_) return false;
  if (Qdvdv.rows() != dimv_) return false;
  if (Qdvdv.cols() != dimv_) return false;
  if (Qxu.rows() != 2*dimv_) return false;
  if (Qxu.cols() != dimu_) return false;
  if (Quu.rows() != dimu_) return false;
  if (Quu.cols() != dimu_) return false;
  if (fx.size() != 2*dimv_) return false;
  if (hx.size() != 2*dimv_) return false;
  if (ha.size() != dimv_) return false;
  if (hu.size() != dimu_) return false;
  return true;
}


bool SplitKKTMatrix::isApprox(const SplitKKTMatrix& other) const {
  if (!Fxx.isApprox(other.Fxx)) return false;
  if (!Fvu.isApprox(other.Fvu)) return false;
  if (!Qxx.isApprox(other.Qxx)) return false;
  if (!Qaa.isApprox(other.Qaa)) return false;
  if (!Qdvdv.isApprox(other.Qdvdv)) return false;
  if (!Qxu.isApprox(other.Qxu)) return false;
  if (!Quu.isApprox(other.Quu)) return false;
  if (!Qff().isApprox(other.Qff())) return false;
  if (!Qqf().isApprox(other.Qqf())) return false;
  if (!fx.isApprox(other.fx)) return false;
  Eigen::VectorXd vec(2), other_vec(2);
  vec << Qtt, Qtt_prev;
  other_vec << other.Qtt, other.Qtt_prev;
  if (!vec.isApprox(other_vec)) return false;
  if (!hx.isApprox(other.hx)) return false;
  if (!hu.isApprox(other.hu)) return false;
  if (!ha.isApprox(other.ha)) return false;
  if (!hf().isApprox(other.hf())) return false;
  return true;
}


bool SplitKKTMatrix::hasNaN() const {
  if (Fxx.hasNaN()) return true;
  if (Fvu.hasNaN()) return true;
  if (Qxx.hasNaN()) return true;
  if (Qaa.hasNaN()) return true;
  if (Qdvdv.hasNaN()) return true;
  if (Qxu.hasNaN()) return true;
  if (Quu.hasNaN()) return true;
  if (Qff().hasNaN()) return true;
  if (Qqf().hasNaN()) return true;
  if (fx.hasNaN()) return true;
  Eigen::VectorXd vec(2);
  vec << Qtt, Qtt_prev;
  if (vec.hasNaN()) return true;
  if (hx.hasNaN()) return true;
  if (hu.hasNaN()) return true;
  if (ha.hasNaN()) return true;
  if (hf().hasNaN()) return true;
  return false;
}


void SplitKKTMatrix::setRandom() {
  Fxx.setRandom();
  Fvu.setRandom();
  const Eigen::MatrixXd Qxxuu_seed = Eigen::MatrixXd::Random(dimx_+dimu_, dimx_+dimu_);
  const Eigen::MatrixXd Qxxuu = Qxxuu_seed * Qxxuu_seed.transpose();
  Qxx = Qxxuu.topLeftCorner(dimx_, dimx_);
  Qxu = Qxxuu.topRightCorner(dimx_, dimu_);
  Quu = Qxxuu.bottomRightCorner(dimu_, dimu_);
  const Eigen::MatrixXd Qaaff_seed = Eigen::MatrixXd::Random(dimv_+dimf_, dimv_+dimf_);
  const Eigen::MatrixXd Qaaff = Qaaff_seed * Qaaff_seed.transpose();
  Qaa = Qaaff.topLeftCorner(dimv_, dimv_);
  Qdvdv = Qaa;
  Qff() = Qaaff.bottomRightCorner(dimf_, dimf_);
  Qqf().setRandom();
  fx.setRandom();
  Qtt = Eigen::VectorXd::Random(1)[0];
  Qtt_prev = Eigen::VectorXd::Random(1)[0];
  hx.setRandom();
  hu.setRandom();
  ha.setRandom();
  hf().setRandom();
}


void SplitKKTMatrix::setRandom(const ContactStatus& contact_status) {
  setContactStatus(contact_status);
  setRandom();
}


SplitKKTMatrix SplitKKTMatrix::Random(const Robot& robot) {
  SplitKKTMatrix kkt_matrix(robot);
  kkt_matrix.setRandom();
  return kkt_matrix;
}


SplitKKTMatrix SplitKKTMatrix::Random(const Robot& robot, 
                                      const ContactStatus& contact_status) {
  SplitKKTMatrix kkt_matrix(robot);
  kkt_matrix.setRandom(contact_status);
  return kkt_matrix;
}

void SplitKKTMatrix::disp(std::ostream& os) const {
  os << "SplitKKTMatrix:" << std::endl;
  os << "  Fxx = " << Fxx << std::endl;
  os << "  Fvu = " << Fvu << std::endl;
  os << "  Qxx = " << Qxx << std::endl;
  os << "  Qxu = " << Qxu << std::endl;
  os << "  Quu = " << Quu << std::endl;
  os << "  Qaa = " << Qaa << std::endl;
  if (dimf_ > 0) {
    os << "  Qff = " << Qff() << std::endl;
    os << "  Qqf = " << Qqf() << std::endl;
  }
  os << "  fq = " << fq().transpose() << std::endl;
  os << "  fv = " << fv().transpose() << std::endl;
  os << "  Qtt = " << Qtt << std::endl;
  os << "  Qtt_prev = " << Qtt_prev << std::endl;
  os << "  hq = " << hq().transpose() << std::endl;
  os << "  hv = " << hv().transpose() << std::endl;
  os << "  hu = " << hu.transpose() << std::endl;
  os << "  ha = " << ha.transpose() << std::flush;
  if (dimf_ > 0) {
    os << std::endl;
    os << "  hf = " << hf().transpose() << std::flush;
  }
}


std::ostream& operator<<(std::ostream& os, const SplitKKTMatrix& kkt_matrix) {
  kkt_matrix.disp(os);
  return os;
}

} // namespace robotoc 