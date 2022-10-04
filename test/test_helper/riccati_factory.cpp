#include "riccati_factory.hpp"

#include "robotoc/ocp/time_discretization.hpp"


namespace robotoc {
namespace testhelper {

SplitRiccatiFactorization CreateSplitRiccatiFactorization(const Robot& robot) {
  auto riccati_factorization = SplitRiccatiFactorization::Random(robot);
  const int dimx = 2*robot.dimv();
  Eigen::MatrixXd seed = Eigen::MatrixXd::Random(dimx, dimx);
  riccati_factorization.P = seed * seed.transpose();
  riccati_factorization.xi = 1000.0 * std::abs(riccati_factorization.xi); // scaling to avoid ill-conditioning
  riccati_factorization.rho = std::abs(riccati_factorization.rho);
  return riccati_factorization;
}

} // namespace testhelper
} // namespace robotoc