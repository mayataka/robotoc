#ifndef IDOCP_IMPULSE_RICCATI_FACTORIZER_HPP_
#define IDOCP_IMPULSE_RICCATI_FACTORIZER_HPP_

#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/ocp/riccati_solution.hpp"
#include "idocp/impulse/impulse_kkt_matrix.hpp"
#include "idocp/impulse/impulse_kkt_residual.hpp"


namespace idocp {

class ImpulseRiccatiFactorizer {
public:
  // Constructor.
  // Argments:
  //    robot: The robot model that has been already initialized.
  ImpulseRiccatiFactorizer(const Robot& robot);

  // Default constructor.
  ImpulseRiccatiFactorizer();

  // Destructor.
  ~ImpulseRiccatiFactorizer();
 
  // Use default copy constructor.
  ImpulseRiccatiFactorizer(const ImpulseRiccatiFactorizer&) = default;

  // Use default copy operator.
  ImpulseRiccatiFactorizer& operator=(const ImpulseRiccatiFactorizer&) = default;

  // Use default move constructor.
  ImpulseRiccatiFactorizer(ImpulseRiccatiFactorizer&&) noexcept = default;

  // Use default move assign operator.
  ImpulseRiccatiFactorizer& operator=(ImpulseRiccatiFactorizer&&) noexcept = default;

  void factorize(const ImpulseKKTMatrix& kkt_matrix,
                 const ImpulseKKTResidual& kkt_residual,
                 const RiccatiSolution& riccati_next, RiccatiSolution& riccati);

  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

private:
  bool has_floating_base_;
  int dimv_;
  static constexpr int kDimFloatingBase = 6;
  Eigen::MatrixXd AtPqq_, AtPqv_, AtPvq_, AtPvv_;

};

} // namespace idocp

#include "idocp/impulse/impulse_riccati_factorizer.hxx"

#endif // IDOCP_IMPULSE_RICCATI_FACTORIZER_HPP_ 