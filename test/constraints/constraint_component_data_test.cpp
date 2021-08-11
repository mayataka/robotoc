#include <gtest/gtest.h>

#include "idocp/constraints/constraint_component_data.hpp"

namespace idocp {

class ConstraintComponentDataTest : public ::testing::Test {
protected:
  virtual void SetUp() {
  }

  virtual void TearDown() {
  }

};


TEST_F(ConstraintComponentDataTest, constructor) {
  const int dimc = 5;
  const double barrier = 0.01;
  ConstraintComponentData data(dimc, barrier);
  EXPECT_EQ(data.slack.size(), dimc);
  EXPECT_EQ(data.dual.size(), dimc);
  EXPECT_EQ(data.residual.size(), dimc);
  EXPECT_EQ(data.cmpl.size(), dimc);
  EXPECT_EQ(data.dslack.size(), dimc);
  EXPECT_EQ(data.ddual.size(), dimc);
  EXPECT_DOUBLE_EQ(data.log_barrier, 0.0);
  EXPECT_EQ(data.dimc(), dimc);
}


TEST_F(ConstraintComponentDataTest, err) {
  const int dimc = 5;
  const double barrier = 0.01;
  ConstraintComponentData data(dimc, barrier);
  data.residual.setRandom();
  data.cmpl.setRandom();

  const double err = data.KKTError();
  const double err_ref = data.residual.squaredNorm() + data.cmpl.squaredNorm();
  EXPECT_DOUBLE_EQ(err, err_ref);

  const double vio = data.constraintViolation();
  const double vio_ref = data.residual.template lpNorm<1>();
  EXPECT_DOUBLE_EQ(vio, vio_ref);
}

} // namespace idocp


int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}