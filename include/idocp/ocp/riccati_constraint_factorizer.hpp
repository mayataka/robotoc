#ifndef IDOCP_RICCATI_CONSTRAINT_FACTORIZER_HPP_
#define IDOCP_RICCATI_CONSTRAINT_FACTORIZER_HPP_

#include "Eigen/Core"
#include "Eigen/LU"

#include "idocp/robot/robot.hpp"
#include "idocp/ocp/kkt_matrix.hpp"
#include "idocp/ocp/kkt_residual.hpp"
#include "idocp/ocp/split_direction.hpp"
#include "idocp/ocp/riccati_solution.hpp"
#include "idocp/ocp/riccati_constraint_solution.hpp"


namespace idocp {

class RiccatiConstraintFactorizer {
public:
  // Constructor.
  // Argments:
  //    robot: The robot model that has been already initialized.
  RiccatiConstraintFactorizer(const Robot& robot, const int N, const int nproc);

  // Default constructor.
  RiccatiConstraintFactorizer();

  // Destructor.
  ~RiccatiConstraintFactorizer();
 
  // Default copy constructor.
  RiccatiConstraintFactorizer(const RiccatiConstraintFactorizer&) = default;

  // Default copy operator.
  RiccatiConstraintFactorizer& operator=(const RiccatiConstraintFactorizer&) = default;

  // Default move constructor.
  RiccatiConstraintFactorizer(RiccatiConstraintFactorizer&&) noexcept = default;

  // Default move assign operator.
  RiccatiConstraintFactorizer& operator=(RiccatiConstraintFactorizer&&) noexcept = default;

  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

  void solve(
      const std::vector<RiccatiSolution>& riccati_impulse, 
      const std::vector<RiccatiConstraintSolution>& riccati_constraint_impulse);

  template <typename MatrixType1, typename MatrixType2, typename MatrixType3>
  static void computeENEtinv(Eigen::LLT<Eigen::MatrixXd>& llt,
                             const Eigen::MatrixBase<MatrixType1>& Eq,
                             const Eigen::MatrixBase<MatrixType2>& Nqq,
                             const Eigen::MatrixBase<MatrixType3>& ENEtinv);

  void invertBlockLowerTriangular(
      const std::vector<RiccatiSolution>& riccati_impulse, 
      const std::vector<RiccatiConstraintSolution>& riccati_constraint_impulse);

  void invertBlockUpperTriangular();
  void invertBlockDiagonal();


private:
  std::vector<Eigen::LLT<Eigen::MatrixXd>> llts_;
  std::vector<ConstrainedRiccatiSolution> sol_;
  int N_, nproc_;

};

} // namespace idocp

#include "idocp/ocp/riccati_constraint_factorizer.hxx"

#endif // IDOCP_RICCATI_CONSTRAINT_FACTORIZER_HPP_ 