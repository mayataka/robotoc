#include <string>

#include <gtest/gtest.h>
#include "Eigen/Core"
#include "Eigen/LU"

#include "idocp/robot/robot.hpp"
#include "idocp/ocp/riccati_matrix_inverter.hpp"


namespace idocp {

class RiccatiMatrixInverterTest : public ::testing::Test {
protected:
  virtual void SetUp() {
    srand((unsigned int) time(0));
    std::random_device rnd;
    fixed_base_urdf_ = "../urdf/iiwa14/iiwa14.urdf";
    floating_base_urdf_ = "../urdf/anymal/anymal.urdf";
    fixed_base_robot_ = Robot(fixed_base_urdf_);
    floating_base_robot_ = Robot(floating_base_urdf_);
    dtau_ = std::abs(Eigen::VectorXd::Random(1)[0]);
  }

  virtual void TearDown() {
  }

  double dtau_;
  std::string fixed_base_urdf_, floating_base_urdf_;
  Robot fixed_base_robot_, floating_base_robot_;
};


TEST_F(RiccatiMatrixInverterTest, fixed_base_without_contacts) {
  const int dimv = fixed_base_robot_.dimv();
  Eigen::MatrixXd Qqa = Eigen::MatrixXd::Random(dimv, dimv);
  Eigen::MatrixXd Qva = Eigen::MatrixXd::Random(dimv, dimv);
  Eigen::MatrixXd Qaa = Eigen::MatrixXd::Random(dimv, dimv);
  Eigen::VectorXd la = Eigen::VectorXd::Random(dimv);
  Eigen::MatrixXd Kaq = Eigen::MatrixXd::Zero(dimv, dimv);
  Eigen::MatrixXd Kav = Eigen::MatrixXd::Zero(dimv, dimv);
  Eigen::VectorXd ka = Eigen::VectorXd::Zero(dimv);
  Qaa.triangularView<Eigen::StrictlyLower>() 
      = Qaa.transpose().triangularView<Eigen::StrictlyLower>();
  // Makes Qaa semi positive define
  Qaa = Qaa * Qaa.transpose();
  // Adds identity matrix to make Qaa sufficiently positive define
  Qaa.noalias() += Eigen::MatrixXd::Identity(dimv, dimv);
  while (Qaa.determinant() == 0) {
    Qaa = Eigen::MatrixXd::Random(dimv, dimv);
    Qaa.triangularView<Eigen::StrictlyLower>() 
        = Qaa.transpose().triangularView<Eigen::StrictlyLower>();
    Qaa = Qaa * Qaa.transpose();
    Qaa.noalias() += Eigen::MatrixXd::Identity(dimv, dimv);
  }
  RiccatiMatrixInverter inverter(fixed_base_robot_);
  inverter.invert(Qqa, Qva, Qaa, la, Kaq, Kav, ka);
  const Eigen::MatrixXd Qaa_inv = Qaa.inverse();
  const Eigen::MatrixXd Kaq_ref = - Qaa_inv * Qqa.transpose();
  const Eigen::MatrixXd Kav_ref = - Qaa_inv * Qva.transpose();
  const Eigen::VectorXd ka_ref = - Qaa_inv * la;
  EXPECT_TRUE(Kaq.isApprox(Kaq_ref));
  EXPECT_TRUE(Kav.isApprox(Kav_ref));
  EXPECT_TRUE(ka.isApprox(ka_ref));
  std::cout << "Kaq error:" << std::endl;
  std::cout << Kaq - Kaq_ref << std::endl;
  std::cout << std::endl;
  std::cout << "Kav error:" << std::endl;
  std::cout << Kav - Kav_ref << std::endl;
  std::cout << std::endl;
  std::cout << "ka error:" << std::endl;
  std::cout << ka - ka_ref << std::endl;
  std::cout << std::endl;
}


TEST_F(RiccatiMatrixInverterTest, fixed_base_with_contacts) {
  const int contact_frame = 18;
  const std::vector<int> contact_frames = {contact_frame};
  fixed_base_robot_ = Robot(fixed_base_urdf_, contact_frames, 0, 0);
  const std::vector<bool> active_contacts = {true};
  fixed_base_robot_.setActiveContacts(active_contacts);
  const int dimv = fixed_base_robot_.dimv();
  const int dimf = fixed_base_robot_.max_dimf();
  Eigen::MatrixXd gen_mat = Eigen::MatrixXd::Random(dimv+dimf, dimv+dimf);
  gen_mat.triangularView<Eigen::StrictlyLower>() 
      = gen_mat.transpose().triangularView<Eigen::StrictlyLower>();
  // Makes pos_mat semi positive define
  Eigen::MatrixXd pos_mat = gen_mat * gen_mat.transpose();
  // Adds identity matrix to make pos_mat sufficiently positive define
  pos_mat.noalias() += Eigen::MatrixXd::Identity(dimv+dimf, dimv+dimf);
  while (pos_mat.determinant() == 0) {
    gen_mat = Eigen::MatrixXd::Random(dimv+dimf, dimv+dimf);
    gen_mat.triangularView<Eigen::StrictlyLower>() 
        = gen_mat.transpose().triangularView<Eigen::StrictlyLower>();
    pos_mat = gen_mat * gen_mat.transpose();
    pos_mat.noalias() += Eigen::MatrixXd::Identity(dimv+dimf, dimv+dimf);
  }
  Eigen::MatrixXd Qqa = Eigen::MatrixXd::Random(dimv, dimv);
  Eigen::MatrixXd Qva = Eigen::MatrixXd::Random(dimv, dimv);
  Eigen::MatrixXd Qaa = pos_mat.block(0, 0, dimv, dimv);
  Eigen::MatrixXd Qqf = Eigen::MatrixXd::Random(dimv, dimf);
  Eigen::MatrixXd Qvf = Eigen::MatrixXd::Random(dimv, dimf);
  Eigen::MatrixXd Qaf = pos_mat.block(0, dimv, dimv, dimf);
  Eigen::MatrixXd Qff = pos_mat.block(dimv, dimv, dimf, dimf);
  Eigen::MatrixXd Cq = Eigen::MatrixXd::Random(dimf, dimv);
  Eigen::MatrixXd Cv = Eigen::MatrixXd::Random(dimf, dimv);
  Eigen::MatrixXd Ca = Eigen::MatrixXd::Random(dimf, dimv);
  Eigen::VectorXd la = Eigen::VectorXd::Random(dimv);
  Eigen::VectorXd lf = Eigen::VectorXd::Random(dimf);
  Eigen::VectorXd C_res = Eigen::VectorXd::Random(dimf);
  Eigen::MatrixXd Kaq = Eigen::MatrixXd::Zero(dimv, dimv);
  Eigen::MatrixXd Kav = Eigen::MatrixXd::Zero(dimv, dimv);
  Eigen::MatrixXd Kfq = Eigen::MatrixXd::Zero(dimf, dimv);
  Eigen::MatrixXd Kfv = Eigen::MatrixXd::Zero(dimf, dimv);
  Eigen::MatrixXd Kmuq = Eigen::MatrixXd::Zero(dimf, dimv);
  Eigen::MatrixXd Kmuv = Eigen::MatrixXd::Zero(dimf, dimv);
  Eigen::VectorXd ka = Eigen::VectorXd::Zero(dimv);
  Eigen::VectorXd kf = Eigen::VectorXd::Zero(dimf);
  Eigen::VectorXd kmu = Eigen::VectorXd::Zero(dimf);
  Eigen::MatrixXd M = Eigen::MatrixXd::Zero(dimv+2*dimf, dimv+2*dimf);
  M.block(0, 0, dimv+dimf, dimv+dimf) = pos_mat;
  M.block(0, dimv+dimf, dimv, dimf) = Ca.transpose();
  M.block(dimv+dimf, 0, dimf, dimv) = Ca;
  M.triangularView<Eigen::StrictlyLower>() 
      = M.transpose().triangularView<Eigen::StrictlyLower>();
  const Eigen::MatrixXd Minv = M.inverse();

  RiccatiMatrixInverter inverter(fixed_base_robot_);
  inverter.setContactStatus(fixed_base_robot_);
  inverter.precompute(Qaf, Qff);
  inverter.invert(Qqa, Qva, Qaa, Qqf, Qvf, Cq, Cv, Ca, la, lf, C_res, Kaq, Kav, 
                  Kfq, Kfv, Kmuq, Kmuv, ka, kf, kmu);

  const Eigen::MatrixXd Kaq_ref = - Minv.block(0, 0, dimv, dimv) * Qqa.transpose()
                                  - Minv.block(0, dimv, dimv, dimf) * Qqf.transpose()
                                  - Minv.block(0, dimv+dimf, dimv, dimf) * Cq;
  const Eigen::MatrixXd Kav_ref = - Minv.block(0, 0, dimv, dimv) * Qva.transpose()
                                  - Minv.block(0, dimv, dimv, dimf) * Qvf.transpose()
                                  - Minv.block(0, dimv+dimf, dimv, dimf) * Cv;
  const Eigen::MatrixXd Kfq_ref = - Minv.block(dimv, 0, dimf, dimv) * Qqa.transpose()
                                  - Minv.block(dimv, dimv, dimf, dimf) * Qqf.transpose()
                                  - Minv.block(dimv, dimv+dimf, dimf, dimf) * Cq;
  const Eigen::MatrixXd Kfv_ref = - Minv.block(dimv, 0, dimf, dimv) * Qva.transpose()
                                  - Minv.block(dimv, dimv, dimf, dimf) * Qvf.transpose()
                                  - Minv.block(dimv, dimv+dimf, dimf, dimf) * Cv;
  const Eigen::MatrixXd Kmuq_ref = - Minv.block(dimv+dimf, 0, dimf, dimv) * Qqa.transpose()
                                   - Minv.block(dimv+dimf, dimv, dimf, dimf) * Qqf.transpose()
                                   - Minv.block(dimv+dimf, dimv+dimf, dimf, dimf) * Cq;
  const Eigen::MatrixXd Kmuv_ref = - Minv.block(dimv+dimf, 0, dimf, dimv) * Qva.transpose()
                                   - Minv.block(dimv+dimf, dimv, dimf, dimf) * Qvf.transpose()
                                   - Minv.block(dimv+dimf, dimv+dimf, dimf, dimf) * Cv;
  const Eigen::VectorXd ka_ref = - Minv.block(0, 0, dimv, dimv) * la
                                 - Minv.block(0, dimv, dimv, dimf) * lf
                                 - Minv.block(0, dimv+dimf, dimv, dimf) * C_res;
  const Eigen::VectorXd kf_ref = - Minv.block(dimv, 0, dimf, dimv) * la
                                 - Minv.block(dimv, dimv, dimf, dimf) * lf
                                 - Minv.block(dimv, dimv+dimf, dimf, dimf) * C_res;
  const Eigen::VectorXd kmu_ref = - Minv.block(dimv+dimf, 0, dimf, dimv) * la
                                  - Minv.block(dimv+dimf, dimv, dimf, dimf) * lf
                                  - Minv.block(dimv+dimf, dimv+dimf, dimf, dimf) * C_res;
  EXPECT_TRUE(Kaq.isApprox(Kaq_ref));
  EXPECT_TRUE(Kav.isApprox(Kav_ref));
  EXPECT_TRUE(Kfq.isApprox(Kfq_ref));
  EXPECT_TRUE(Kfv.isApprox(Kfv_ref));
  EXPECT_TRUE(Kmuq.isApprox(Kmuq_ref));
  EXPECT_TRUE(Kmuv.isApprox(Kmuv_ref));
  EXPECT_TRUE(ka.isApprox(ka_ref));
  EXPECT_TRUE(kf.isApprox(kf_ref));
  EXPECT_TRUE(kmu.isApprox(kmu_ref));
  std::cout << "Kaq error:" << std::endl;
  std::cout << Kaq - Kaq_ref << std::endl;
  std::cout << std::endl;
  std::cout << "Kav error:" << std::endl;
  std::cout << Kav - Kav_ref << std::endl;
  std::cout << std::endl;
  std::cout << "Kfq error:" << std::endl;
  std::cout << Kfq - Kfq_ref << std::endl;
  std::cout << std::endl;
  std::cout << "Kfv error:" << std::endl;
  std::cout << Kfv - Kfv_ref << std::endl;
  std::cout << std::endl;
  std::cout << "Kmuq error:" << std::endl;
  std::cout << Kmuq - Kmuq_ref << std::endl;
  std::cout << std::endl;
  std::cout << "Kmuv error:" << std::endl;
  std::cout << Kmuv - Kmuv_ref << std::endl;
  std::cout << std::endl;
  std::cout << "ka error:" << std::endl;
  std::cout << ka - ka_ref << std::endl;
  std::cout << std::endl;
  std::cout << "kf error:" << std::endl;
  std::cout << kf - kf_ref << std::endl;
  std::cout << std::endl;
  std::cout << "kmu error:" << std::endl;
  std::cout << kmu - kmu_ref << std::endl;
  std::cout << std::endl;
}


TEST_F(RiccatiMatrixInverterTest, floating_base_without_contacts) {
  const int dimv = floating_base_robot_.dimv();
  const int dimf = floating_base_robot_.max_dimf();
  const int dimc = dimf + floating_base_robot_.dim_passive();
  ASSERT_TRUE(dimf == 0);
  ASSERT_TRUE(dimc == 6);
  Eigen::MatrixXd gen_mat = Eigen::MatrixXd::Random(dimv, dimv);
  gen_mat.triangularView<Eigen::StrictlyLower>() 
      = gen_mat.transpose().triangularView<Eigen::StrictlyLower>();
  // Makes pos_mat semi positive define
  Eigen::MatrixXd pos_mat = gen_mat * gen_mat.transpose();
  // Adds identity matrix to make pos_mat sufficiently positive define
  pos_mat.noalias() += Eigen::MatrixXd::Identity(dimv, dimv);
  while (pos_mat.determinant() == 0) {
    gen_mat = Eigen::MatrixXd::Random(dimv, dimv);
    gen_mat.triangularView<Eigen::StrictlyLower>() 
        = gen_mat.transpose().triangularView<Eigen::StrictlyLower>();
    pos_mat = gen_mat * gen_mat.transpose();
    pos_mat.noalias() += Eigen::MatrixXd::Identity(dimv, dimv);
  }
  Eigen::MatrixXd Qqa = Eigen::MatrixXd::Random(dimv, dimv);
  Eigen::MatrixXd Qva = Eigen::MatrixXd::Random(dimv, dimv);
  Eigen::MatrixXd Qaa = pos_mat.block(0, 0, dimv, dimv);
  Eigen::MatrixXd Cq = Eigen::MatrixXd::Random(dimc, dimv);
  Eigen::MatrixXd Cv = Eigen::MatrixXd::Random(dimc, dimv);
  Eigen::MatrixXd Ca = Eigen::MatrixXd::Random(dimc, dimv);
  Eigen::VectorXd la = Eigen::VectorXd::Random(dimv);
  Eigen::VectorXd C_res = Eigen::VectorXd::Random(dimc);
  Eigen::MatrixXd Kaq = Eigen::MatrixXd::Zero(dimv, dimv);
  Eigen::MatrixXd Kav = Eigen::MatrixXd::Zero(dimv, dimv);
  Eigen::MatrixXd Kmuq = Eigen::MatrixXd::Zero(dimc, dimv);
  Eigen::MatrixXd Kmuv = Eigen::MatrixXd::Zero(dimc, dimv);
  Eigen::VectorXd ka = Eigen::VectorXd::Zero(dimv);
  Eigen::VectorXd kmu = Eigen::VectorXd::Zero(dimc);
  Eigen::MatrixXd M = Eigen::MatrixXd::Zero(dimv+dimc, dimv+dimc);
  M.block(0, 0, dimv, dimv) = pos_mat;
  M.block(0, dimv, dimv, dimc) = Ca.transpose();
  M.block(dimv, 0, dimc, dimv) = Ca;
  M.triangularView<Eigen::StrictlyLower>() 
      = M.transpose().triangularView<Eigen::StrictlyLower>();
  const Eigen::MatrixXd Minv = M.inverse();

  RiccatiMatrixInverter inverter(floating_base_robot_);
  inverter.setContactStatus(floating_base_robot_);
  inverter.invert(Qqa, Qva, Qaa, Cq, Cv, Ca, la, C_res, Kaq, Kav, Kmuq, Kmuv, 
                  ka, kmu);

  const Eigen::MatrixXd Kaq_ref = - Minv.block(0, 0, dimv, dimv) * Qqa.transpose()
                                  - Minv.block(0, dimv, dimv, dimc) * Cq;
  const Eigen::MatrixXd Kav_ref = - Minv.block(0, 0, dimv, dimv) * Qva.transpose()
                                  - Minv.block(0, dimv, dimv, dimc) * Cv;
  const Eigen::MatrixXd Kmuq_ref = - Minv.block(dimv, 0, dimc, dimv) * Qqa.transpose()
                                   - Minv.block(dimv, dimv, dimc, dimc) * Cq;
  const Eigen::MatrixXd Kmuv_ref = - Minv.block(dimv, 0, dimc, dimv) * Qva.transpose()
                                   - Minv.block(dimv, dimv, dimc, dimc) * Cv;
  const Eigen::VectorXd ka_ref = - Minv.block(0, 0, dimv, dimv) * la
                                 - Minv.block(0, dimv, dimv, dimc) * C_res;
  const Eigen::VectorXd kmu_ref = - Minv.block(dimv, 0, dimc, dimv) * la
                                  - Minv.block(dimv, dimv, dimc, dimc) * C_res;
  EXPECT_TRUE(Kaq.isApprox(Kaq_ref));
  EXPECT_TRUE(Kav.isApprox(Kav_ref));
  EXPECT_TRUE(Kmuq.isApprox(Kmuq_ref));
  EXPECT_TRUE(Kmuv.isApprox(Kmuv_ref));
  EXPECT_TRUE(ka.isApprox(ka_ref));
  EXPECT_TRUE(kmu.isApprox(kmu_ref));
  std::cout << "Kaq error:" << std::endl;
  std::cout << Kaq - Kaq_ref << std::endl;
  std::cout << std::endl;
  std::cout << "Kav error:" << std::endl;
  std::cout << Kav - Kav_ref << std::endl;
  std::cout << std::endl;
  std::cout << "Kmuq error:" << std::endl;
  std::cout << Kmuq - Kmuq_ref << std::endl;
  std::cout << std::endl;
  std::cout << "Kmuv error:" << std::endl;
  std::cout << Kmuv - Kmuv_ref << std::endl;
  std::cout << std::endl;
  std::cout << "ka error:" << std::endl;
  std::cout << ka - ka_ref << std::endl;
  std::cout << std::endl;
  std::cout << "kmu error:" << std::endl;
  std::cout << kmu - kmu_ref << std::endl;
  std::cout << std::endl;
}


TEST_F(RiccatiMatrixInverterTest, floating_base_with_contacts) {
  const std::vector<int> contact_frames = {14, 24, 34, 44};
  floating_base_robot_ = Robot(floating_base_urdf_, contact_frames, 0, 0);
  const std::vector<bool> active_contacts = {true, true, true, true};
  floating_base_robot_.setActiveContacts(active_contacts);
  const int dimv = floating_base_robot_.dimv();
  const int dimf = floating_base_robot_.max_dimf();
  const int dimc = dimf + floating_base_robot_.dim_passive();
  const int dim_passive = floating_base_robot_.dim_passive();
  Eigen::MatrixXd gen_mat = Eigen::MatrixXd::Random(dimv+dimf, dimv+dimf);
  gen_mat.triangularView<Eigen::StrictlyLower>() 
      = gen_mat.transpose().triangularView<Eigen::StrictlyLower>();
  // Makes pos_mat semi positive define
  Eigen::MatrixXd pos_mat = gen_mat * gen_mat.transpose();
  // Adds identity matrix to make pos_mat sufficiently positive define
  pos_mat.noalias() += Eigen::MatrixXd::Identity(dimv+dimf, dimv+dimf);
  while (pos_mat.determinant() == 0) {
    gen_mat = Eigen::MatrixXd::Random(dimv+dimf, dimv+dimf);
    gen_mat.triangularView<Eigen::StrictlyLower>() 
        = gen_mat.transpose().triangularView<Eigen::StrictlyLower>();
    pos_mat = gen_mat * gen_mat.transpose();
    pos_mat.noalias() += Eigen::MatrixXd::Identity(dimv+dimf, dimv+dimf);
  }
  Eigen::MatrixXd Qqa = Eigen::MatrixXd::Random(dimv, dimv);
  Eigen::MatrixXd Qva = Eigen::MatrixXd::Random(dimv, dimv);
  Eigen::MatrixXd Qaa = pos_mat.block(0, 0, dimv, dimv);
  Eigen::MatrixXd Qqf = Eigen::MatrixXd::Random(dimv, dimf);
  Eigen::MatrixXd Qvf = Eigen::MatrixXd::Random(dimv, dimf);
  Eigen::MatrixXd Qaf = pos_mat.block(0, dimv, dimv, dimf);
  Eigen::MatrixXd Qff = pos_mat.block(dimv, dimv, dimf, dimf);
  Eigen::MatrixXd Cq = Eigen::MatrixXd::Random(dimc, dimv);
  Eigen::MatrixXd Cv = Eigen::MatrixXd::Random(dimc, dimv);
  Eigen::MatrixXd Ca = Eigen::MatrixXd::Random(dimc, dimv);
  Eigen::MatrixXd Cf = Eigen::MatrixXd::Random(dim_passive, dimf);
  Eigen::VectorXd la = Eigen::VectorXd::Random(dimv);
  Eigen::VectorXd lf = Eigen::VectorXd::Random(dimf);
  Eigen::VectorXd C_res = Eigen::VectorXd::Random(dimc);
  Eigen::MatrixXd Kaq = Eigen::MatrixXd::Zero(dimv, dimv);
  Eigen::MatrixXd Kav = Eigen::MatrixXd::Zero(dimv, dimv);
  Eigen::MatrixXd Kfq = Eigen::MatrixXd::Zero(dimf, dimv);
  Eigen::MatrixXd Kfv = Eigen::MatrixXd::Zero(dimf, dimv);
  Eigen::MatrixXd Kmuq = Eigen::MatrixXd::Zero(dimc, dimv);
  Eigen::MatrixXd Kmuv = Eigen::MatrixXd::Zero(dimc, dimv);
  Eigen::VectorXd ka = Eigen::VectorXd::Zero(dimv);
  Eigen::VectorXd kf = Eigen::VectorXd::Zero(dimf);
  Eigen::VectorXd kmu = Eigen::VectorXd::Zero(dimc);
  Eigen::MatrixXd M = Eigen::MatrixXd::Zero(dimv+dimf+dimc, dimv+dimf+dimc);
  M.block(0, 0, dimv+dimf, dimv+dimf) = pos_mat;
  M.block(0, dimv+dimf, dimv, dimc) = Ca.transpose();
  M.block(dimv, dimv+dimf, dimf, dim_passive) = Cf.transpose();
  M.triangularView<Eigen::StrictlyLower>() 
      = M.transpose().triangularView<Eigen::StrictlyLower>();
  const Eigen::MatrixXd Minv = M.inverse();

  RiccatiMatrixInverter inverter(floating_base_robot_);
  inverter.setContactStatus(floating_base_robot_);
  inverter.precompute(Qaf, Qff);
  inverter.invert(Qqa, Qva, Qaa, Qqf, Qvf, Cq, Cv, Ca, Cf, la, lf, C_res, 
                  Kaq, Kav, Kfq, Kfv, Kmuq, Kmuv, ka, kf, kmu);
  const Eigen::MatrixXd Kaq_ref = - Minv.block(0, 0, dimv, dimv) * Qqa.transpose()
                                  - Minv.block(0, dimv, dimv, dimf) * Qqf.transpose()
                                  - Minv.block(0, dimv+dimf, dimv, dimc) * Cq;
  const Eigen::MatrixXd Kav_ref = - Minv.block(0, 0, dimv, dimv) * Qva.transpose()
                                  - Minv.block(0, dimv, dimv, dimf) * Qvf.transpose()
                                  - Minv.block(0, dimv+dimf, dimv, dimc) * Cv;
  const Eigen::MatrixXd Kfq_ref = - Minv.block(dimv, 0, dimf, dimv) * Qqa.transpose()
                                  - Minv.block(dimv, dimv, dimf, dimf) * Qqf.transpose()
                                  - Minv.block(dimv, dimv+dimf, dimf, dimc) * Cq;
  const Eigen::MatrixXd Kfv_ref = - Minv.block(dimv, 0, dimf, dimv) * Qva.transpose()
                                  - Minv.block(dimv, dimv, dimf, dimf) * Qvf.transpose()
                                  - Minv.block(dimv, dimv+dimf, dimf, dimc) * Cv;
  const Eigen::MatrixXd Kmuq_ref = - Minv.block(dimv+dimf, 0, dimc, dimv) * Qqa.transpose()
                                   - Minv.block(dimv+dimf, dimv, dimc, dimf) * Qqf.transpose()
                                   - Minv.block(dimv+dimf, dimv+dimf, dimc, dimc) * Cq;
  const Eigen::MatrixXd Kmuv_ref = - Minv.block(dimv+dimf, 0, dimc, dimv) * Qva.transpose()
                                   - Minv.block(dimv+dimf, dimv, dimc, dimf) * Qvf.transpose()
                                   - Minv.block(dimv+dimf, dimv+dimf, dimc, dimc) * Cv;
  const Eigen::VectorXd ka_ref = - Minv.block(0, 0, dimv, dimv) * la
                                 - Minv.block(0, dimv, dimv, dimf) * lf
                                 - Minv.block(0, dimv+dimf, dimv, dimc) * C_res;
  const Eigen::VectorXd kf_ref = - Minv.block(dimv, 0, dimf, dimv) * la
                                 - Minv.block(dimv, dimv, dimf, dimf) * lf
                                 - Minv.block(dimv, dimv+dimf, dimf, dimc) * C_res;
  const Eigen::VectorXd kmu_ref = - Minv.block(dimv+dimf, 0, dimc, dimv) * la
                                  - Minv.block(dimv+dimf, dimv, dimc, dimf) * lf
                                  - Minv.block(dimv+dimf, dimv+dimf, dimc, dimc) * C_res;
  EXPECT_TRUE(Kaq.isApprox(Kaq_ref));
  EXPECT_TRUE(Kav.isApprox(Kav_ref));
  EXPECT_TRUE(Kfq.isApprox(Kfq_ref));
  EXPECT_TRUE(Kfv.isApprox(Kfv_ref));
  EXPECT_TRUE(Kmuq.isApprox(Kmuq_ref));
  EXPECT_TRUE(Kmuv.isApprox(Kmuv_ref));
  EXPECT_TRUE(ka.isApprox(ka_ref));
  EXPECT_TRUE(kf.isApprox(kf_ref));
  EXPECT_TRUE(kmu.isApprox(kmu_ref));
  std::cout << "Kaq error:" << std::endl;
  std::cout << Kaq - Kaq_ref << std::endl;
  std::cout << std::endl;
  std::cout << "Kav error:" << std::endl;
  std::cout << Kav - Kav_ref << std::endl;
  std::cout << std::endl;
  std::cout << "Kfq error:" << std::endl;
  std::cout << Kfq - Kfq_ref << std::endl;
  std::cout << std::endl;
  std::cout << "Kfv error:" << std::endl;
  std::cout << Kfv - Kfv_ref << std::endl;
  std::cout << std::endl;
  std::cout << "Kmuq error:" << std::endl;
  std::cout << Kmuq - Kmuq_ref << std::endl;
  std::cout << std::endl;
  std::cout << "Kmuv error:" << std::endl;
  std::cout << Kmuv - Kmuv_ref << std::endl;
  std::cout << std::endl;
  std::cout << "ka error:" << std::endl;
  std::cout << ka - ka_ref << std::endl;
  std::cout << std::endl;
  std::cout << "kf error:" << std::endl;
  std::cout << kf - kf_ref << std::endl;
  std::cout << std::endl;
  std::cout << "kmu error:" << std::endl;
  std::cout << kmu - kmu_ref << std::endl;
  std::cout << std::endl;
}

} // namespace idocp


int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}