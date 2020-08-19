#ifndef IDOCP_CHOLESKY_HPP_
#define IDOCP_CHOLESKY_HPP_

#include <assert.h>

#include "Eigen/Core"
#include "Eigen/LU"
#include "Eigen/Cholesky"

#include "idocp/ocp/kkt_matrix.hpp"


namespace idocp {

class Cholesky {
public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

  Cholesky(const Robot& robot) 
    : max_dimKKT_(5*robot.dimv()+2*robot.max_dimf()+robot.dim_passive()),
      L(Eigen::Matrix::Zero(max_dimKKT_, max_dimKKT_), 
      D(Eigen::Matrix::Zero(max_dimKKT_, max_dimKKT_) {
  }

  Cholesky() 
    : kkt_composition_(),
      kkt_matrix_() {
  }

  ~Cholesky() {
  }

  Cholesky(const Cholesky&) = default;

  Cholesky& operator=(const Cholesky&) = default;
 
  Cholesky(Cholesky&&) noexcept = default;

  Cholesky& operator=(Cholesky&&) noexcept = default;

  void invert(const KKTMatrix& kkt_matrix);

private:
  Eigen::MatrixXd L, D;
  int max_dimKKT_;

};

} // namespace idocp 


#endif // IDOCP_CHOLESKY_HPP_