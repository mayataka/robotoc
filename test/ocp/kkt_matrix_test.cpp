#include <string>

#include <gtest/gtest.h>

#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/ocp/kkt_matrix.hpp"


namespace idocp {

class KKTMatrixTest : public ::testing::Test {
protected:
  virtual void SetUp() {
    srand((unsigned int) time(0));
    std::random_device rnd;
    fixed_base_urdf_ = "../urdf/iiwa14/iiwa14.urdf";
    floating_base_urdf_ = "../urdf/anymal/anymal.urdf";
  }

  virtual void TearDown() {
  }

  double dtau_;
  std::string fixed_base_urdf_, floating_base_urdf_;
};


TEST_F(KKTMatrixTest, fixed_base) {
  std::vector<int> contact_frames = {18};
  const std::vector<double> mu = {std::abs(Eigen::VectorXd::Random(1)[0])};
  Robot robot(fixed_base_urdf_, contact_frames, mu, 0, 0);
  std::random_device rnd;
  std::vector<bool> contact_status = {rnd()%2==0};
  robot.setContactStatus(contact_status);
  KKTMatrix matrix(robot);
  matrix.setContactStatus(robot);
  const int dimv = robot.dimv();
  const int dimf = robot.dimf();
  const int dim_passive = robot.dim_passive();
  const int dimc = robot.dim_passive() + robot.dimf();
  EXPECT_EQ(matrix.dimf(), dimf);
  EXPECT_EQ(matrix.dimc(), dimc);
  EXPECT_EQ(matrix.dimKKT(), 5*dimv+dimf+dimc);
  const Eigen::MatrixXd Ca = Eigen::MatrixXd::Random(dimc, dimv);
  const Eigen::MatrixXd Cf = Eigen::MatrixXd::Random(dimc, dimf);
  const Eigen::MatrixXd Cq = Eigen::MatrixXd::Random(dimc, dimv);
  const Eigen::MatrixXd Cv = Eigen::MatrixXd::Random(dimc, dimv);
  const Eigen::MatrixXd Qaa = Eigen::MatrixXd::Random(dimv, dimv);
  const Eigen::MatrixXd Qaf = Eigen::MatrixXd::Random(dimv, dimf);
  const Eigen::MatrixXd Qaq = Eigen::MatrixXd::Random(dimv, dimv);
  const Eigen::MatrixXd Qav = Eigen::MatrixXd::Random(dimv, dimv);
  const Eigen::MatrixXd Qfa = Eigen::MatrixXd::Random(dimf, dimv);
  const Eigen::MatrixXd Qff = Eigen::MatrixXd::Random(dimf, dimf);
  const Eigen::MatrixXd Qfq = Eigen::MatrixXd::Random(dimf, dimv);
  const Eigen::MatrixXd Qfv = Eigen::MatrixXd::Random(dimf, dimv);
  const Eigen::MatrixXd Qqa = Eigen::MatrixXd::Random(dimv, dimv);
  const Eigen::MatrixXd Qqf = Eigen::MatrixXd::Random(dimv, dimf);
  const Eigen::MatrixXd Qqq = Eigen::MatrixXd::Random(dimv, dimv);
  const Eigen::MatrixXd Qqv = Eigen::MatrixXd::Random(dimv, dimv);
  const Eigen::MatrixXd Qva = Eigen::MatrixXd::Random(dimv, dimv);
  const Eigen::MatrixXd Qvf = Eigen::MatrixXd::Random(dimv, dimf);
  const Eigen::MatrixXd Qvq = Eigen::MatrixXd::Random(dimv, dimv);
  const Eigen::MatrixXd Qvv = Eigen::MatrixXd::Random(dimv, dimv);
  matrix.Ca() = Ca;
  matrix.Cf() = Cf;
  matrix.Cq() = Cq;
  matrix.Cv() = Cv;
  matrix.Qaa() = Qaa;
  matrix.Qaf() = Qaf;
  matrix.Qaq() = Qaq;
  matrix.Qav() = Qav;
  matrix.Qfa() = Qfa;
  matrix.Qff() = Qff;
  matrix.Qfq() = Qfq;
  matrix.Qfv() = Qfv;
  matrix.Qqa() = Qqa;
  matrix.Qqf() = Qqf;
  matrix.Qqq() = Qqq;
  matrix.Qqv() = Qqv;
  matrix.Qva() = Qva;
  matrix.Qvf() = Qvf;
  matrix.Qvq() = Qvq;
  matrix.Qvv() = Qvv;
  EXPECT_TRUE(matrix.constraintsJacobian().block(0,           0, dimc, dimv).isApprox(Ca));
  EXPECT_TRUE(matrix.constraintsJacobian().block(0,        dimv, dimc, dimf).isApprox(Cf));
  EXPECT_TRUE(matrix.constraintsJacobian().block(0,   dimv+dimf, dimc, dimv).isApprox(Cq));
  EXPECT_TRUE(matrix.constraintsJacobian().block(0, 2*dimv+dimf, dimc, dimv).isApprox(Cv));
  EXPECT_TRUE(matrix.Ca().topRows(robot.dim_passive()).isApprox(matrix.Ca_floating_base()));
  EXPECT_TRUE(matrix.Cf().topRows(robot.dim_passive()).isApprox(matrix.Cf_floating_base()));
  EXPECT_TRUE(matrix.Cq().topRows(robot.dim_passive()).isApprox(matrix.Cq_floating_base()));
  EXPECT_TRUE(matrix.Cv().topRows(robot.dim_passive()).isApprox(matrix.Cv_floating_base()));
  EXPECT_TRUE(matrix.Ca().bottomRows(robot.dimf()).isApprox(matrix.Ca_contacts()));
  EXPECT_TRUE(matrix.Cf().bottomRows(robot.dimf()).isApprox(matrix.Cf_contacts()));
  EXPECT_TRUE(matrix.Cq().bottomRows(robot.dimf()).isApprox(matrix.Cq_contacts()));
  EXPECT_TRUE(matrix.Cv().bottomRows(robot.dimf()).isApprox(matrix.Cv_contacts()));
  EXPECT_TRUE(matrix.costHessian().block(          0,           0, dimv, dimv).isApprox(Qaa));
  EXPECT_TRUE(matrix.costHessian().block(          0,        dimv, dimv, dimf).isApprox(Qaf));
  EXPECT_TRUE(matrix.costHessian().block(          0,   dimv+dimf, dimv, dimv).isApprox(Qaq));
  EXPECT_TRUE(matrix.costHessian().block(          0, 2*dimv+dimf, dimv, dimv).isApprox(Qav));
  EXPECT_TRUE(matrix.costHessian().block(       dimv,           0, dimf, dimv).isApprox(Qfa));
  EXPECT_TRUE(matrix.costHessian().block(       dimv,        dimv, dimf, dimf).isApprox(Qff));
  EXPECT_TRUE(matrix.costHessian().block(       dimv,   dimv+dimf, dimf, dimv).isApprox(Qfq));
  EXPECT_TRUE(matrix.costHessian().block(       dimv, 2*dimv+dimf, dimf, dimv).isApprox(Qfv));
  EXPECT_TRUE(matrix.costHessian().block(  dimv+dimf,           0, dimv, dimv).isApprox(Qqa));
  EXPECT_TRUE(matrix.costHessian().block(  dimv+dimf,        dimv, dimv, dimf).isApprox(Qqf));
  EXPECT_TRUE(matrix.costHessian().block(  dimv+dimf,   dimv+dimf, dimv, dimv).isApprox(Qqq));
  EXPECT_TRUE(matrix.costHessian().block(  dimv+dimf, 2*dimv+dimf, dimv, dimv).isApprox(Qqv));
  EXPECT_TRUE(matrix.costHessian().block(2*dimv+dimf,           0, dimv, dimv).isApprox(Qva));
  EXPECT_TRUE(matrix.costHessian().block(2*dimv+dimf,        dimv, dimv, dimf).isApprox(Qvf));
  EXPECT_TRUE(matrix.costHessian().block(2*dimv+dimf,   dimv+dimf, dimv, dimv).isApprox(Qvq));
  EXPECT_TRUE(matrix.costHessian().block(2*dimv+dimf, 2*dimv+dimf, dimv, dimv).isApprox(Qvv));
  EXPECT_EQ(matrix.Qxx().rows(), 2*dimv);
  EXPECT_EQ(matrix.Qxx().cols(), 2*dimv);
  EXPECT_TRUE(matrix.Qxx().block(   0,    0, dimv, dimv).isApprox(Qqq));
  EXPECT_TRUE(matrix.Qxx().block(   0, dimv, dimv, dimv).isApprox(Qqv));
  EXPECT_TRUE(matrix.Qxx().block(dimv,    0, dimv, dimv).isApprox(Qvq));
  EXPECT_TRUE(matrix.Qxx().block(dimv, dimv, dimv, dimv).isApprox(Qvv));
  EXPECT_TRUE(matrix.Qxx().block(   0,    0, dimv, dimv).isApprox(Qqq));
  EXPECT_TRUE(matrix.Qxx().block(   0, dimv, dimv, dimv).isApprox(Qqv));
  EXPECT_TRUE(matrix.Qxx().block(dimv,    0, dimv, dimv).isApprox(Qvq));
  EXPECT_TRUE(matrix.Qxx().block(dimv, dimv, dimv, dimv).isApprox(Qvv));
  EXPECT_EQ(matrix.Qafaf().rows(), dimv+dimf);
  EXPECT_EQ(matrix.Qafaf().cols(), dimv+dimf);
  EXPECT_TRUE(matrix.Qafaf().block(   0,    0, dimv, dimv).isApprox(Qaa));
  EXPECT_TRUE(matrix.Qafaf().block(   0, dimv, dimv, dimf).isApprox(Qaf));
  EXPECT_TRUE(matrix.Qafaf().block(dimv,    0, dimf, dimv).isApprox(Qfa));
  EXPECT_TRUE(matrix.Qafaf().block(dimv, dimv, dimf, dimf).isApprox(Qff));
  EXPECT_EQ(matrix.Qafqv().rows(), dimv+dimf);
  EXPECT_EQ(matrix.Qafqv().cols(), 2*dimv);
  EXPECT_TRUE(matrix.Qafqv().block(   0,    0, dimv, dimv).isApprox(Qaq));
  EXPECT_TRUE(matrix.Qafqv().block(   0, dimv, dimv, dimv).isApprox(Qav));
  EXPECT_TRUE(matrix.Qafqv().block(dimv,    0, dimf, dimv).isApprox(Qfq));
  EXPECT_TRUE(matrix.Qafqv().block(dimv, dimv, dimf, dimv).isApprox(Qfv));
  EXPECT_EQ(matrix.Cqv().rows(), dimc);
  EXPECT_EQ(matrix.Cqv().cols(), 2*dimv);
  EXPECT_TRUE(matrix.Cqv().block(0, 0, dimc, dimv).isApprox(Cq));
  EXPECT_TRUE(matrix.Cqv().block(0, dimv, dimc, dimv).isApprox(Cv));
  EXPECT_EQ(matrix.Caf().rows(), dimc);
  EXPECT_EQ(matrix.Caf().cols(), dimv+dimf);
  EXPECT_TRUE(matrix.Caf().block(0, 0, dimc, dimv).isApprox(Ca));
  EXPECT_TRUE(matrix.Caf().block(0, dimv, dimc, dimf).isApprox(Cf));
  matrix.setZero();
  EXPECT_TRUE(matrix.costHessian().isZero());
  EXPECT_TRUE(matrix.constraintsJacobian().isZero());
}


TEST_F(KKTMatrixTest, invert_fixed_base) {
  std::vector<int> contact_frames = {18};
  const std::vector<double> mu = {std::abs(Eigen::VectorXd::Random(1)[0])};
  Robot robot(fixed_base_urdf_, contact_frames, mu, 0, 0);
  std::random_device rnd;
  std::vector<bool> contact_status = {rnd()%2==0};
  robot.setContactStatus(contact_status);
  KKTMatrix matrix(robot);
  matrix.setContactStatus(robot);
  const int dimv = robot.dimv();
  const int dimx = 2*robot.dimv();
  const int dimf = robot.dimf();
  const int dim_passive = robot.dim_passive();
  const int dimc = robot.dim_passive() + robot.dimf();
  const int dimQ = 3*robot.dimv() + robot.dimf();
  const Eigen::MatrixXd Q_seed_mat = Eigen::MatrixXd::Random(dimQ, dimQ);
  const Eigen::MatrixXd Q_mat = Q_seed_mat * Q_seed_mat.transpose() + Eigen::MatrixXd::Identity(dimQ, dimQ);
  const Eigen::MatrixXd Jc_mat = Eigen::MatrixXd::Random(dimc, dimQ);
  matrix.costHessian() = Q_mat;
  matrix.constraintsJacobian() = Jc_mat;
  const double dtau = std::abs(Eigen::VectorXd::Random(1)[0]);
  const int dimKKT = 5*robot.dimv()+robot.dim_passive()+2*robot.dimf();
  Eigen::MatrixXd kkt_mat_ref = Eigen::MatrixXd::Zero(dimKKT, dimKKT);
  kkt_mat_ref.bottomRightCorner(dimQ, dimQ) = Q_mat;
  kkt_mat_ref.block(dimx, dimx+dimc, dimc, dimQ) = Jc_mat;
  kkt_mat_ref.block(0, dimx+dimc+dimv+dimf, dimv, dimv) 
      = -1 * Eigen::MatrixXd::Identity(dimv, dimv);
  kkt_mat_ref.block(0, dimx+dimc+2*dimv+dimf, dimv, dimv) 
      = dtau * Eigen::MatrixXd::Identity(dimv, dimv);
  kkt_mat_ref.block(dimv, dimx+dimc, dimv, dimv) 
      = dtau * Eigen::MatrixXd::Identity(dimv, dimv);
  kkt_mat_ref.block(dimv, dimx+dimc+2*dimv+dimf, dimv, dimv) 
      = -1 * Eigen::MatrixXd::Identity(dimv, dimv);
  kkt_mat_ref.triangularView<Eigen::StrictlyLower>() 
      = kkt_mat_ref.transpose().triangularView<Eigen::StrictlyLower>();
  std::cout << kkt_mat_ref << std::endl;
  const Eigen::MatrixXd kkt_mat_inv_ref = kkt_mat_ref.inverse();
  Eigen::MatrixXd kkt_mat_inv = Eigen::MatrixXd::Zero(dimKKT, dimKKT);
  matrix.invert(dtau, kkt_mat_inv);
  EXPECT_TRUE(kkt_mat_inv.isApprox(kkt_mat_inv_ref, 1.0e-08));
  std::cout << "error l2 norm = " << (kkt_mat_inv - kkt_mat_inv_ref).lpNorm<2>() << std::endl;
}


TEST_F(KKTMatrixTest, floating_base) {
  std::vector<int> contact_frames = {14, 24, 34, 44};
  std::vector<double> mu;
  for (int i=0; i<contact_frames.size(); ++i) {
    mu.push_back(std::abs(Eigen::VectorXd::Random(1)[0]));
  }
  Robot robot(floating_base_urdf_, contact_frames, mu, 0, 0);
  std::random_device rnd;
  std::vector<bool> contact_status;
  for (const auto frame : contact_frames) {
    contact_status.push_back(rnd()%2==0);
  }
  robot.setContactStatus(contact_status);
  KKTMatrix matrix(robot);
  matrix.setContactStatus(robot);
  const int dimv = robot.dimv();
  const int dimf = robot.dimf();
  const int dim_passive = robot.dim_passive();
  const int dimc = robot.dim_passive() + robot.dimf();
  EXPECT_EQ(matrix.dimf(), dimf);
  EXPECT_EQ(matrix.dimc(), dimc);
  EXPECT_EQ(matrix.dimKKT(), 5*dimv+dimf+dimc);
  const Eigen::MatrixXd Ca = Eigen::MatrixXd::Random(dimc, dimv);
  const Eigen::MatrixXd Cf = Eigen::MatrixXd::Random(dimc, dimf);
  const Eigen::MatrixXd Cq = Eigen::MatrixXd::Random(dimc, dimv);
  const Eigen::MatrixXd Cv = Eigen::MatrixXd::Random(dimc, dimv);
  const Eigen::MatrixXd Qaa = Eigen::MatrixXd::Random(dimv, dimv);
  const Eigen::MatrixXd Qaf = Eigen::MatrixXd::Random(dimv, dimf);
  const Eigen::MatrixXd Qaq = Eigen::MatrixXd::Random(dimv, dimv);
  const Eigen::MatrixXd Qav = Eigen::MatrixXd::Random(dimv, dimv);
  const Eigen::MatrixXd Qfa = Eigen::MatrixXd::Random(dimf, dimv);
  const Eigen::MatrixXd Qff = Eigen::MatrixXd::Random(dimf, dimf);
  const Eigen::MatrixXd Qfq = Eigen::MatrixXd::Random(dimf, dimv);
  const Eigen::MatrixXd Qfv = Eigen::MatrixXd::Random(dimf, dimv);
  const Eigen::MatrixXd Qqa = Eigen::MatrixXd::Random(dimv, dimv);
  const Eigen::MatrixXd Qqf = Eigen::MatrixXd::Random(dimv, dimf);
  const Eigen::MatrixXd Qqq = Eigen::MatrixXd::Random(dimv, dimv);
  const Eigen::MatrixXd Qqv = Eigen::MatrixXd::Random(dimv, dimv);
  const Eigen::MatrixXd Qva = Eigen::MatrixXd::Random(dimv, dimv);
  const Eigen::MatrixXd Qvf = Eigen::MatrixXd::Random(dimv, dimf);
  const Eigen::MatrixXd Qvq = Eigen::MatrixXd::Random(dimv, dimv);
  const Eigen::MatrixXd Qvv = Eigen::MatrixXd::Random(dimv, dimv);
  matrix.Ca() = Ca;
  matrix.Cf() = Cf;
  matrix.Cq() = Cq;
  matrix.Cv() = Cv;
  matrix.Qaa() = Qaa;
  matrix.Qaf() = Qaf;
  matrix.Qaq() = Qaq;
  matrix.Qav() = Qav;
  matrix.Qfa() = Qfa;
  matrix.Qff() = Qff;
  matrix.Qfq() = Qfq;
  matrix.Qfv() = Qfv;
  matrix.Qqa() = Qqa;
  matrix.Qqf() = Qqf;
  matrix.Qqq() = Qqq;
  matrix.Qqv() = Qqv;
  matrix.Qva() = Qva;
  matrix.Qvf() = Qvf;
  matrix.Qvq() = Qvq;
  matrix.Qvv() = Qvv;
  EXPECT_TRUE(matrix.constraintsJacobian().block(0,           0, dimc, dimv).isApprox(Ca));
  EXPECT_TRUE(matrix.constraintsJacobian().block(0,        dimv, dimc, dimf).isApprox(Cf));
  EXPECT_TRUE(matrix.constraintsJacobian().block(0,   dimv+dimf, dimc, dimv).isApprox(Cq));
  EXPECT_TRUE(matrix.constraintsJacobian().block(0, 2*dimv+dimf, dimc, dimv).isApprox(Cv));
  EXPECT_TRUE(matrix.Ca().topRows(robot.dim_passive()).isApprox(matrix.Ca_floating_base()));
  EXPECT_TRUE(matrix.Cf().topRows(robot.dim_passive()).isApprox(matrix.Cf_floating_base()));
  EXPECT_TRUE(matrix.Cq().topRows(robot.dim_passive()).isApprox(matrix.Cq_floating_base()));
  EXPECT_TRUE(matrix.Cv().topRows(robot.dim_passive()).isApprox(matrix.Cv_floating_base()));
  EXPECT_TRUE(matrix.Ca().bottomRows(robot.dimf()).isApprox(matrix.Ca_contacts()));
  EXPECT_TRUE(matrix.Cf().bottomRows(robot.dimf()).isApprox(matrix.Cf_contacts()));
  EXPECT_TRUE(matrix.Cq().bottomRows(robot.dimf()).isApprox(matrix.Cq_contacts()));
  EXPECT_TRUE(matrix.Cv().bottomRows(robot.dimf()).isApprox(matrix.Cv_contacts()));
  EXPECT_TRUE(matrix.costHessian().block(          0,           0, dimv, dimv).isApprox(Qaa));
  EXPECT_TRUE(matrix.costHessian().block(          0,        dimv, dimv, dimf).isApprox(Qaf));
  EXPECT_TRUE(matrix.costHessian().block(          0,   dimv+dimf, dimv, dimv).isApprox(Qaq));
  EXPECT_TRUE(matrix.costHessian().block(          0, 2*dimv+dimf, dimv, dimv).isApprox(Qav));
  EXPECT_TRUE(matrix.costHessian().block(       dimv,           0, dimf, dimv).isApprox(Qfa));
  EXPECT_TRUE(matrix.costHessian().block(       dimv,        dimv, dimf, dimf).isApprox(Qff));
  EXPECT_TRUE(matrix.costHessian().block(       dimv,   dimv+dimf, dimf, dimv).isApprox(Qfq));
  EXPECT_TRUE(matrix.costHessian().block(       dimv, 2*dimv+dimf, dimf, dimv).isApprox(Qfv));
  EXPECT_TRUE(matrix.costHessian().block(  dimv+dimf,           0, dimv, dimv).isApprox(Qqa));
  EXPECT_TRUE(matrix.costHessian().block(  dimv+dimf,        dimv, dimv, dimf).isApprox(Qqf));
  EXPECT_TRUE(matrix.costHessian().block(  dimv+dimf,   dimv+dimf, dimv, dimv).isApprox(Qqq));
  EXPECT_TRUE(matrix.costHessian().block(  dimv+dimf, 2*dimv+dimf, dimv, dimv).isApprox(Qqv));
  EXPECT_TRUE(matrix.costHessian().block(2*dimv+dimf,           0, dimv, dimv).isApprox(Qva));
  EXPECT_TRUE(matrix.costHessian().block(2*dimv+dimf,        dimv, dimv, dimf).isApprox(Qvf));
  EXPECT_TRUE(matrix.costHessian().block(2*dimv+dimf,   dimv+dimf, dimv, dimv).isApprox(Qvq));
  EXPECT_TRUE(matrix.costHessian().block(2*dimv+dimf, 2*dimv+dimf, dimv, dimv).isApprox(Qvv));
  EXPECT_EQ(matrix.Qxx().rows(), 2*dimv);
  EXPECT_EQ(matrix.Qxx().cols(), 2*dimv);
  EXPECT_TRUE(matrix.Qxx().block(   0,    0, dimv, dimv).isApprox(Qqq));
  EXPECT_TRUE(matrix.Qxx().block(   0, dimv, dimv, dimv).isApprox(Qqv));
  EXPECT_TRUE(matrix.Qxx().block(dimv,    0, dimv, dimv).isApprox(Qvq));
  EXPECT_TRUE(matrix.Qxx().block(dimv, dimv, dimv, dimv).isApprox(Qvv));
  EXPECT_TRUE(matrix.Qxx().block(   0,    0, dimv, dimv).isApprox(Qqq));
  EXPECT_TRUE(matrix.Qxx().block(   0, dimv, dimv, dimv).isApprox(Qqv));
  EXPECT_TRUE(matrix.Qxx().block(dimv,    0, dimv, dimv).isApprox(Qvq));
  EXPECT_TRUE(matrix.Qxx().block(dimv, dimv, dimv, dimv).isApprox(Qvv));
  EXPECT_EQ(matrix.Qafaf().rows(), dimv+dimf);
  EXPECT_EQ(matrix.Qafaf().cols(), dimv+dimf);
  EXPECT_TRUE(matrix.Qafaf().block(   0,    0, dimv, dimv).isApprox(Qaa));
  EXPECT_TRUE(matrix.Qafaf().block(   0, dimv, dimv, dimf).isApprox(Qaf));
  EXPECT_TRUE(matrix.Qafaf().block(dimv,    0, dimf, dimv).isApprox(Qfa));
  EXPECT_TRUE(matrix.Qafaf().block(dimv, dimv, dimf, dimf).isApprox(Qff));
  EXPECT_EQ(matrix.Qafqv().rows(), dimv+dimf);
  EXPECT_EQ(matrix.Qafqv().cols(), 2*dimv);
  EXPECT_TRUE(matrix.Qafqv().block(   0,    0, dimv, dimv).isApprox(Qaq));
  EXPECT_TRUE(matrix.Qafqv().block(   0, dimv, dimv, dimv).isApprox(Qav));
  EXPECT_TRUE(matrix.Qafqv().block(dimv,    0, dimf, dimv).isApprox(Qfq));
  EXPECT_TRUE(matrix.Qafqv().block(dimv, dimv, dimf, dimv).isApprox(Qfv));
  EXPECT_EQ(matrix.Cqv().rows(), dimc);
  EXPECT_EQ(matrix.Cqv().cols(), 2*dimv);
  EXPECT_TRUE(matrix.Cqv().block(0, 0, dimc, dimv).isApprox(Cq));
  EXPECT_TRUE(matrix.Cqv().block(0, dimv, dimc, dimv).isApprox(Cv));
  EXPECT_EQ(matrix.Caf().rows(), dimc);
  EXPECT_EQ(matrix.Caf().cols(), dimv+dimf);
  EXPECT_TRUE(matrix.Caf().block(0, 0, dimc, dimv).isApprox(Ca));
  EXPECT_TRUE(matrix.Caf().block(0, dimv, dimc, dimf).isApprox(Cf));
  matrix.setZero();
  EXPECT_TRUE(matrix.costHessian().isZero());
  EXPECT_TRUE(matrix.constraintsJacobian().isZero());
}


TEST_F(KKTMatrixTest, invert_floating_base) {
  std::vector<int> contact_frames = {14, 24, 34, 44};
  std::vector<double> mu;
  for (int i=0; i<contact_frames.size(); ++i) {
    mu.push_back(std::abs(Eigen::VectorXd::Random(1)[0]));
  }
  Robot robot(floating_base_urdf_, contact_frames, mu, 0, 0);
  std::random_device rnd;
  std::vector<bool> contact_status;
  for (const auto frame : contact_frames) {
    contact_status.push_back(rnd()%2==0);
  }
  robot.setContactStatus(contact_status);
  KKTMatrix matrix(robot);
  matrix.setContactStatus(robot);
  const int dimv = robot.dimv();
  const int dimx = 2*robot.dimv();
  const int dimf = robot.dimf();
  const int dim_passive = robot.dim_passive();
  const int dimc = robot.dim_passive() + robot.dimf();
  const int dimQ = 3*robot.dimv() + robot.dimf();
  const Eigen::MatrixXd Q_seed_mat = Eigen::MatrixXd::Random(dimQ, dimQ);
  const Eigen::MatrixXd Q_mat = Q_seed_mat * Q_seed_mat.transpose() + Eigen::MatrixXd::Identity(dimQ, dimQ);
  const Eigen::MatrixXd Jc_mat = Eigen::MatrixXd::Random(dimc, dimQ);
  matrix.costHessian() = Q_mat;
  matrix.constraintsJacobian() = Jc_mat;
  const Eigen::MatrixXd Fqq_seed_mat = Eigen::MatrixXd::Random(6, 6);
  const Eigen::MatrixXd Fqq_mat = Fqq_seed_mat * Fqq_seed_mat.transpose();
  matrix.Fqq = -1 * Eigen::MatrixXd::Identity(dimv, dimv);
  matrix.Fqq.topLeftCorner(6, 6) = - Fqq_mat;
  const double dtau = std::abs(Eigen::VectorXd::Random(1)[0]);
  const int dimKKT = 5*robot.dimv()+robot.dim_passive()+2*robot.dimf();
  Eigen::MatrixXd kkt_mat_ref = Eigen::MatrixXd::Zero(dimKKT, dimKKT);
  kkt_mat_ref.bottomRightCorner(dimQ, dimQ) = Q_mat;
  kkt_mat_ref.block(dimx, dimx+dimc, dimc, dimQ) = Jc_mat;
  kkt_mat_ref.block(0, dimx+dimc+dimv+dimf, dimv, dimv) 
      = -1 * Eigen::MatrixXd::Identity(dimv, dimv);
  kkt_mat_ref.block(0, dimx+dimc+dimv+dimf, dimv, dimv).topLeftCorner(6, 6)
      = - Fqq_mat;
  kkt_mat_ref.block(0, dimx+dimc+2*dimv+dimf, dimv, dimv) 
      = dtau * Eigen::MatrixXd::Identity(dimv, dimv);
  kkt_mat_ref.block(dimv, dimx+dimc, dimv, dimv) 
      = dtau * Eigen::MatrixXd::Identity(dimv, dimv);
  kkt_mat_ref.block(dimv, dimx+dimc+2*dimv+dimf, dimv, dimv) 
      = -1 * Eigen::MatrixXd::Identity(dimv, dimv);
  kkt_mat_ref.triangularView<Eigen::StrictlyLower>() 
      = kkt_mat_ref.transpose().triangularView<Eigen::StrictlyLower>();
  std::cout << kkt_mat_ref << std::endl;
  const Eigen::MatrixXd kkt_mat_inv_ref = kkt_mat_ref.inverse();
  Eigen::MatrixXd kkt_mat_inv = Eigen::MatrixXd::Zero(dimKKT, dimKKT);
  matrix.invert(dtau, kkt_mat_inv);
  EXPECT_TRUE(kkt_mat_inv.isApprox(kkt_mat_inv_ref, 1.0e-08));
  std::cout << "error l2 norm = " << (kkt_mat_inv - kkt_mat_inv_ref).lpNorm<2>() << std::endl;
}

} // namespace idocp


int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}