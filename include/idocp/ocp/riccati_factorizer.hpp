#ifndef IDOCP_RICCATI_FACTORIZER_HPP_
#define IDOCP_RICCATI_FACTORIZER_HPP_

#include "Eigen/Core"
#include "Eigen/LU"

#include "idocp/robot/robot.hpp"
#include "idocp/ocp/kkt_matrix.hpp"
#include "idocp/ocp/kkt_residual.hpp"
#include "idocp/ocp/split_direction.hpp"
#include "idocp/ocp/riccati_factorization.hpp"


namespace idocp {

class RiccatiFactorizer {
public:
  // Constructor.
  // Argments:
  //    robot: The robot model that has been already initialized.
  RiccatiFactorizer(const Robot& robot);

  // Default constructor.
  RiccatiFactorizer();

  // Destructor.
  ~RiccatiFactorizer();
 
  // Use default copy constructor.
  RiccatiFactorizer(const RiccatiFactorizer&) = default;

  // Use default copy operator.
  RiccatiFactorizer& operator=(const RiccatiFactorizer&) = default;

  // Use default move constructor.
  RiccatiFactorizer(RiccatiFactorizer&&) noexcept = default;

  // Use default move assign operator.
  RiccatiFactorizer& operator=(RiccatiFactorizer&&) noexcept = default;

  void backwardRiccatiRecursion(const RiccatiFactorization& riccati_next, 
                                const double dtau, KKTMatrix& kkt_matrix, 
                                KKTResidual& kkt_residual, 
                                RiccatiFactorization& riccati);

  void forwardRiccatiRecursion(const KKTMatrix& kkt_matrix, 
                               const KKTResidual& kkt_residual,
                               const SplitDirection& d, const double dtau,
                               SplitDirection& d_next) const;

  void factorizeStateConstraintParallel(const KKTMatrix& kkt_matrix, 
                                        const KKTResidual& kkt_residual, 
                                        const double dtau,
                                        RiccatiFactorization& riccati);

  static void factorizeStateConstraintSerial(RiccatiFactorization& riccati);

  void factorizeStateConstraintSerial(const RiccatiFactorization& riccati,
                                      RiccatiFactorization& riccati_next);

  static void computeCostateDirection(const RiccatiFactorization& riccati, 
                                      SplitDirection& d);

  static void computeControlInputDirection(const RiccatiFactorization& riccati, 
                                           SplitDirection& d);

  void factorizeKKTMatrix(const RiccatiFactorization& riccati_next, 
                          const double dtau, KKTMatrix& kkt_matrix, 
                          KKTResidual& kkt_residual);

  void computeFeedbackGainAndFeedforward(const KKTMatrix& kkt_matrix, 
                                         const KKTResidual& kkt_residual,
                                         RiccatiFactorization& riccati);

  void factorizeRiccatiFactorization(const RiccatiFactorization& riccati_next, 
                                     const double dtau, 
                                     const KKTMatrix& kkt_matrix, 
                                     const KKTResidual& kkt_residual,
                                     RiccatiFactorization& riccati);

  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

private:
  bool has_floating_base_;
  int dimv_, dimu_;
  static constexpr int kDimFloatingBase = 6;
  Eigen::LLT<Eigen::MatrixXd> llt_;
  Eigen::MatrixXd AtPqq_, AtPqv_, AtPvq_, AtPvv_, BtPq_, BtPv_, GK_, 
                  GinvBt_, NApBKt_;
};

} // namespace idocp

#include "idocp/ocp/riccati_factorizer.hxx"

#endif // IDOCP_RICCATI_FACTORIZER_HPP_