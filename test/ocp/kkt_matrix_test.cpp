#include <string>

#include <gtest/gtest.h>

#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/ocp/kkt_composition.hpp"
#include "idocp/ocp/kkt_matrix.hpp"


namespace idocp {

class KKTMatrixTest : public ::testing::Test {
protected:
  virtual void SetUp() {
    srand((unsigned int) time(0));
    std::random_device rnd;
    fixed_base_urdf_ = "../urdf/iiwa14/iiwa14.urdf";
    floating_base_urdf_ = "../urdf/iiwa14/iiwa14.urdf";
  }

  virtual void TearDown() {
  }

  double dtau_;
  std::string fixed_base_urdf_, floating_base_urdf_;
};


TEST_F(KKTMatrixTest, fixed_base) {
  std::vector<int> contact_frames = {18};
  Robot robot(fixed_base_urdf_, contact_frames, 0, 0);
  std::random_device rnd;
  std::vector<bool> contact_status = {rnd()%2==0};
  robot.setContactStatus(contact_status);
  KKTComposition composition(robot);
  composition.set(robot);
  KKTMatrix matrix(robot);
  matrix.setContactStatus(robot);
  const Eigen::MatrixXd Fqq = Eigen::MatrixXd::Random(composition.Fq_size(), composition.Fq_size());
  const Eigen::MatrixXd Fqv = Eigen::MatrixXd::Random(composition.Fq_size(), composition.Fv_size());
  const Eigen::MatrixXd Fvq = Eigen::MatrixXd::Random(composition.Fv_size(), composition.Fq_size());
  const Eigen::MatrixXd Fvv = Eigen::MatrixXd::Random(composition.Fv_size(), composition.Fv_size());

  const Eigen::MatrixXd Ca = Eigen::MatrixXd::Random(composition.C_size(), composition.Qa_size());
  const Eigen::MatrixXd Cf = Eigen::MatrixXd::Random(composition.C_size(), composition.Qf_size());
  const Eigen::MatrixXd Cq = Eigen::MatrixXd::Random(composition.C_size(), composition.Qq_size());
  const Eigen::MatrixXd Cv = Eigen::MatrixXd::Random(composition.C_size(), composition.Qv_size());
  const Eigen::MatrixXd Qaa = Eigen::MatrixXd::Random(composition.Qa_size(), composition.Qa_size());
  const Eigen::MatrixXd Qaf = Eigen::MatrixXd::Random(composition.Qa_size(), composition.Qf_size());
  const Eigen::MatrixXd Qaq = Eigen::MatrixXd::Random(composition.Qa_size(), composition.Qq_size());
  const Eigen::MatrixXd Qav = Eigen::MatrixXd::Random(composition.Qa_size(), composition.Qv_size());

  const Eigen::MatrixXd Qfa = Eigen::MatrixXd::Random(composition.Qf_size(), composition.Qa_size());
  const Eigen::MatrixXd Qff = Eigen::MatrixXd::Random(composition.Qf_size(), composition.Qf_size());
  const Eigen::MatrixXd Qfq = Eigen::MatrixXd::Random(composition.Qf_size(), composition.Qq_size());
  const Eigen::MatrixXd Qfv = Eigen::MatrixXd::Random(composition.Qf_size(), composition.Qv_size());

  const Eigen::MatrixXd Qqa = Eigen::MatrixXd::Random(composition.Qq_size(), composition.Qa_size());
  const Eigen::MatrixXd Qqf = Eigen::MatrixXd::Random(composition.Qq_size(), composition.Qf_size());
  const Eigen::MatrixXd Qqq = Eigen::MatrixXd::Random(composition.Qq_size(), composition.Qq_size());
  const Eigen::MatrixXd Qqv = Eigen::MatrixXd::Random(composition.Qq_size(), composition.Qv_size());

  const Eigen::MatrixXd Qva = Eigen::MatrixXd::Random(composition.Qv_size(), composition.Qa_size());
  const Eigen::MatrixXd Qvf = Eigen::MatrixXd::Random(composition.Qv_size(), composition.Qf_size());
  const Eigen::MatrixXd Qvq = Eigen::MatrixXd::Random(composition.Qv_size(), composition.Qq_size());
  const Eigen::MatrixXd Qvv = Eigen::MatrixXd::Random(composition.Qv_size(), composition.Qv_size());
  matrix.Fqq() = Fqq;
  matrix.Fqv() = Fqv;
  matrix.Fvq() = Fvq;
  matrix.Fvv() = Fvv;
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
  EXPECT_TRUE(matrix.KKT_matrix().block(composition.Fq_begin(), composition.Qq_begin(), composition.Fq_size(), composition.Qq_size()).isApprox(Fqq));
  EXPECT_TRUE(matrix.KKT_matrix().block(composition.Fq_begin(), composition.Qv_begin(), composition.Fq_size(), composition.Qv_size()).isApprox(Fqv));
  EXPECT_TRUE(matrix.KKT_matrix().block(composition.Fv_begin(), composition.Qq_begin(), composition.Fv_size(), composition.Qq_size()).isApprox(Fvq));
  EXPECT_TRUE(matrix.KKT_matrix().block(composition.Fv_begin(), composition.Qv_begin(), composition.Fv_size(), composition.Qv_size()).isApprox(Fvv));
  EXPECT_TRUE(matrix.KKT_matrix().block(composition.C_begin(), composition.Qa_begin(), composition.C_size(), composition.Qa_size()).isApprox(Ca));
  EXPECT_TRUE(matrix.KKT_matrix().block(composition.C_begin(), composition.Qf_begin(), composition.C_size(), composition.Qf_size()).isApprox(Cf));
  EXPECT_TRUE(matrix.KKT_matrix().block(composition.C_begin(), composition.Qq_begin(), composition.C_size(), composition.Qq_size()).isApprox(Cq));
  EXPECT_TRUE(matrix.KKT_matrix().block(composition.C_begin(), composition.Qv_begin(), composition.C_size(), composition.Qv_size()).isApprox(Cv));
  EXPECT_TRUE(matrix.KKT_matrix().block(composition.Qa_begin(), composition.Qa_begin(), composition.Qa_size(), composition.Qa_size()).isApprox(Qaa));
  EXPECT_TRUE(matrix.KKT_matrix().block(composition.Qa_begin(), composition.Qf_begin(), composition.Qa_size(), composition.Qf_size()).isApprox(Qaf));
  EXPECT_TRUE(matrix.KKT_matrix().block(composition.Qa_begin(), composition.Qq_begin(), composition.Qa_size(), composition.Qq_size()).isApprox(Qaq));
  EXPECT_TRUE(matrix.KKT_matrix().block(composition.Qa_begin(), composition.Qv_begin(), composition.Qa_size(), composition.Qv_size()).isApprox(Qav));
  EXPECT_TRUE(matrix.KKT_matrix().block(composition.Qf_begin(), composition.Qa_begin(), composition.Qf_size(), composition.Qa_size()).isApprox(Qfa));
  EXPECT_TRUE(matrix.KKT_matrix().block(composition.Qf_begin(), composition.Qf_begin(), composition.Qf_size(), composition.Qf_size()).isApprox(Qff));
  EXPECT_TRUE(matrix.KKT_matrix().block(composition.Qf_begin(), composition.Qq_begin(), composition.Qf_size(), composition.Qq_size()).isApprox(Qfq));
  EXPECT_TRUE(matrix.KKT_matrix().block(composition.Qf_begin(), composition.Qv_begin(), composition.Qf_size(), composition.Qv_size()).isApprox(Qfv));
  EXPECT_TRUE(matrix.KKT_matrix().block(composition.Qq_begin(), composition.Qa_begin(), composition.Qq_size(), composition.Qa_size()).isApprox(Qqa));
  EXPECT_TRUE(matrix.KKT_matrix().block(composition.Qq_begin(), composition.Qf_begin(), composition.Qq_size(), composition.Qf_size()).isApprox(Qqf));
  EXPECT_TRUE(matrix.KKT_matrix().block(composition.Qq_begin(), composition.Qq_begin(), composition.Qq_size(), composition.Qq_size()).isApprox(Qqq));
  EXPECT_TRUE(matrix.KKT_matrix().block(composition.Qq_begin(), composition.Qv_begin(), composition.Qq_size(), composition.Qv_size()).isApprox(Qqv));
  EXPECT_TRUE(matrix.KKT_matrix().block(composition.Qv_begin(), composition.Qa_begin(), composition.Qv_size(), composition.Qa_size()).isApprox(Qva));
  EXPECT_TRUE(matrix.KKT_matrix().block(composition.Qv_begin(), composition.Qf_begin(), composition.Qv_size(), composition.Qf_size()).isApprox(Qvf));
  EXPECT_TRUE(matrix.KKT_matrix().block(composition.Qv_begin(), composition.Qq_begin(), composition.Qv_size(), composition.Qq_size()).isApprox(Qvq));
  EXPECT_TRUE(matrix.KKT_matrix().block(composition.Qv_begin(), composition.Qv_begin(), composition.Qv_size(), composition.Qv_size()).isApprox(Qvv));

  EXPECT_TRUE(matrix.Qxx().block(0, 0, composition.Qq_size(), composition.Qq_size()).isApprox(Qqq));
  EXPECT_TRUE(matrix.Qxx().block(0, composition.Qq_size(), composition.Qq_size(), composition.Qv_size()).isApprox(Qqv));
  EXPECT_TRUE(matrix.Qxx().block(composition.Qq_size(), 0, composition.Qv_size(), composition.Qq_size()).isApprox(Qvq));
  EXPECT_TRUE(matrix.Qxx().block(composition.Qq_size(), composition.Qq_size(), composition.Qv_size(), composition.Qv_size()).isApprox(Qvv));

  matrix.setZero();
  EXPECT_TRUE(matrix.KKT_matrix().isZero());
  EXPECT_EQ(matrix.dimKKT(), composition.dimKKT());
  EXPECT_EQ(matrix.max_dimKKT(), composition.max_dimKKT());
}


TEST_F(KKTMatrixTest, floating_base) {
  std::vector<int> contact_frames = {14, 24, 34, 44};
  Robot robot(fixed_base_urdf_, contact_frames, 0, 0);
  std::random_device rnd;
  std::vector<bool> contact_status;
  for (const auto frame : contact_frames) {
    contact_status.push_back(rnd()%2==0);
  }
  robot.setContactStatus(contact_status);
  KKTComposition composition(robot);
  composition.set(robot);
  KKTMatrix matrix(robot);
  matrix.setContactStatus(robot);
  const Eigen::MatrixXd Fqq = Eigen::MatrixXd::Random(composition.Fq_size(), composition.Fq_size());
  const Eigen::MatrixXd Fqv = Eigen::MatrixXd::Random(composition.Fq_size(), composition.Fv_size());
  const Eigen::MatrixXd Fvq = Eigen::MatrixXd::Random(composition.Fv_size(), composition.Fq_size());
  const Eigen::MatrixXd Fvv = Eigen::MatrixXd::Random(composition.Fv_size(), composition.Fv_size());

  const Eigen::MatrixXd Ca = Eigen::MatrixXd::Random(composition.C_size(), composition.Qa_size());
  const Eigen::MatrixXd Cf = Eigen::MatrixXd::Random(composition.C_size(), composition.Qf_size());
  const Eigen::MatrixXd Cq = Eigen::MatrixXd::Random(composition.C_size(), composition.Qq_size());
  const Eigen::MatrixXd Cv = Eigen::MatrixXd::Random(composition.C_size(), composition.Qv_size());
  const Eigen::MatrixXd Qaa = Eigen::MatrixXd::Random(composition.Qa_size(), composition.Qa_size());
  const Eigen::MatrixXd Qaf = Eigen::MatrixXd::Random(composition.Qa_size(), composition.Qf_size());
  const Eigen::MatrixXd Qaq = Eigen::MatrixXd::Random(composition.Qa_size(), composition.Qq_size());
  const Eigen::MatrixXd Qav = Eigen::MatrixXd::Random(composition.Qa_size(), composition.Qv_size());

  const Eigen::MatrixXd Qfa = Eigen::MatrixXd::Random(composition.Qf_size(), composition.Qa_size());
  const Eigen::MatrixXd Qff = Eigen::MatrixXd::Random(composition.Qf_size(), composition.Qf_size());
  const Eigen::MatrixXd Qfq = Eigen::MatrixXd::Random(composition.Qf_size(), composition.Qq_size());
  const Eigen::MatrixXd Qfv = Eigen::MatrixXd::Random(composition.Qf_size(), composition.Qv_size());

  const Eigen::MatrixXd Qqa = Eigen::MatrixXd::Random(composition.Qq_size(), composition.Qa_size());
  const Eigen::MatrixXd Qqf = Eigen::MatrixXd::Random(composition.Qq_size(), composition.Qf_size());
  const Eigen::MatrixXd Qqq = Eigen::MatrixXd::Random(composition.Qq_size(), composition.Qq_size());
  const Eigen::MatrixXd Qqv = Eigen::MatrixXd::Random(composition.Qq_size(), composition.Qv_size());

  const Eigen::MatrixXd Qva = Eigen::MatrixXd::Random(composition.Qv_size(), composition.Qa_size());
  const Eigen::MatrixXd Qvf = Eigen::MatrixXd::Random(composition.Qv_size(), composition.Qf_size());
  const Eigen::MatrixXd Qvq = Eigen::MatrixXd::Random(composition.Qv_size(), composition.Qq_size());
  const Eigen::MatrixXd Qvv = Eigen::MatrixXd::Random(composition.Qv_size(), composition.Qv_size());
  matrix.Fqq() = Fqq;
  matrix.Fqv() = Fqv;
  matrix.Fvq() = Fvq;
  matrix.Fvv() = Fvv;
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
  EXPECT_TRUE(matrix.KKT_matrix().block(composition.Fq_begin(), composition.Qq_begin(), composition.Fq_size(), composition.Qq_size()).isApprox(Fqq));
  EXPECT_TRUE(matrix.KKT_matrix().block(composition.Fq_begin(), composition.Qv_begin(), composition.Fq_size(), composition.Qv_size()).isApprox(Fqv));
  EXPECT_TRUE(matrix.KKT_matrix().block(composition.Fv_begin(), composition.Qq_begin(), composition.Fv_size(), composition.Qq_size()).isApprox(Fvq));
  EXPECT_TRUE(matrix.KKT_matrix().block(composition.Fv_begin(), composition.Qv_begin(), composition.Fv_size(), composition.Qv_size()).isApprox(Fvv));
  EXPECT_TRUE(matrix.KKT_matrix().block(composition.C_begin(), composition.Qa_begin(), composition.C_size(), composition.Qa_size()).isApprox(Ca));
  EXPECT_TRUE(matrix.KKT_matrix().block(composition.C_begin(), composition.Qf_begin(), composition.C_size(), composition.Qf_size()).isApprox(Cf));
  EXPECT_TRUE(matrix.KKT_matrix().block(composition.C_begin(), composition.Qq_begin(), composition.C_size(), composition.Qq_size()).isApprox(Cq));
  EXPECT_TRUE(matrix.KKT_matrix().block(composition.C_begin(), composition.Qv_begin(), composition.C_size(), composition.Qv_size()).isApprox(Cv));
  EXPECT_TRUE(matrix.KKT_matrix().block(composition.Qa_begin(), composition.Qa_begin(), composition.Qa_size(), composition.Qa_size()).isApprox(Qaa));
  EXPECT_TRUE(matrix.KKT_matrix().block(composition.Qa_begin(), composition.Qf_begin(), composition.Qa_size(), composition.Qf_size()).isApprox(Qaf));
  EXPECT_TRUE(matrix.KKT_matrix().block(composition.Qa_begin(), composition.Qq_begin(), composition.Qa_size(), composition.Qq_size()).isApprox(Qaq));
  EXPECT_TRUE(matrix.KKT_matrix().block(composition.Qa_begin(), composition.Qv_begin(), composition.Qa_size(), composition.Qv_size()).isApprox(Qav));
  EXPECT_TRUE(matrix.KKT_matrix().block(composition.Qf_begin(), composition.Qa_begin(), composition.Qf_size(), composition.Qa_size()).isApprox(Qfa));
  EXPECT_TRUE(matrix.KKT_matrix().block(composition.Qf_begin(), composition.Qf_begin(), composition.Qf_size(), composition.Qf_size()).isApprox(Qff));
  EXPECT_TRUE(matrix.KKT_matrix().block(composition.Qf_begin(), composition.Qq_begin(), composition.Qf_size(), composition.Qq_size()).isApprox(Qfq));
  EXPECT_TRUE(matrix.KKT_matrix().block(composition.Qf_begin(), composition.Qv_begin(), composition.Qf_size(), composition.Qv_size()).isApprox(Qfv));
  EXPECT_TRUE(matrix.KKT_matrix().block(composition.Qq_begin(), composition.Qa_begin(), composition.Qq_size(), composition.Qa_size()).isApprox(Qqa));
  EXPECT_TRUE(matrix.KKT_matrix().block(composition.Qq_begin(), composition.Qf_begin(), composition.Qq_size(), composition.Qf_size()).isApprox(Qqf));
  EXPECT_TRUE(matrix.KKT_matrix().block(composition.Qq_begin(), composition.Qq_begin(), composition.Qq_size(), composition.Qq_size()).isApprox(Qqq));
  EXPECT_TRUE(matrix.KKT_matrix().block(composition.Qq_begin(), composition.Qv_begin(), composition.Qq_size(), composition.Qv_size()).isApprox(Qqv));
  EXPECT_TRUE(matrix.KKT_matrix().block(composition.Qv_begin(), composition.Qa_begin(), composition.Qv_size(), composition.Qa_size()).isApprox(Qva));
  EXPECT_TRUE(matrix.KKT_matrix().block(composition.Qv_begin(), composition.Qf_begin(), composition.Qv_size(), composition.Qf_size()).isApprox(Qvf));
  EXPECT_TRUE(matrix.KKT_matrix().block(composition.Qv_begin(), composition.Qq_begin(), composition.Qv_size(), composition.Qq_size()).isApprox(Qvq));
  EXPECT_TRUE(matrix.KKT_matrix().block(composition.Qv_begin(), composition.Qv_begin(), composition.Qv_size(), composition.Qv_size()).isApprox(Qvv));

  EXPECT_TRUE(matrix.Qxx().block(0, 0, composition.Qq_size(), composition.Qq_size()).isApprox(Qqq));
  EXPECT_TRUE(matrix.Qxx().block(0, composition.Qq_size(), composition.Qq_size(), composition.Qv_size()).isApprox(Qqv));
  EXPECT_TRUE(matrix.Qxx().block(composition.Qq_size(), 0, composition.Qv_size(), composition.Qq_size()).isApprox(Qvq));
  EXPECT_TRUE(matrix.Qxx().block(composition.Qq_size(), composition.Qq_size(), composition.Qv_size(), composition.Qv_size()).isApprox(Qvv));

  matrix.setZero();
  EXPECT_TRUE(matrix.KKT_matrix().isZero());
  EXPECT_EQ(matrix.dimKKT(), composition.dimKKT());
  EXPECT_EQ(matrix.max_dimKKT(), composition.max_dimKKT());
}


} // namespace idocp


int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}