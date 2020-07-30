#ifndef IDOCP_BLOCK_MATRIX_INVERTER_HPP_
#define IDOCP_BLOCK_MATRIX_INVERTER_HPP_

#include <vector>
#include <utility>

#include "Eigen/Core"


namespace idocp {

class RiccatiMatrixInverter {
public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

  RiccatiMatrixInverter(const Robot& robot);

  // Destructor.
  ~RiccatiMatrixInverter();
 
  // Use default copy constructor.
  RiccatiMatrixInverter(const RiccatiMatrixInverter& other) = default;

  // Use default copy operator.
  RiccatiMatrixInverter& operator=(const RiccatiMatrixInverter& other) 
      = default;

  void setContactStatus(const Robot& robot);

  void setPrecomputableMatrices(const Eigen::MatrixXd& Qff, 
                                const Eigen::MatrixXd& Qaf);

  void computeSteteFeedbackGains(const Eigen::MatrixXd& Qaa, 
                                 Eigen::MatrixXd& Kaq, 
                                 Eigen::MatrixXd& Kav);

  void computeSteteFeedbackGains(const Eigen::MatrixXd& Qaa, 
                                 Eigen::MatrixXd& Kaq, 
                                 Eigen::MatrixXd& Kav, 
                                 Eigen::MatrixXd& Kmuq, 
                                 Eigen::MatrixXd& Kmuv);

  void computeSteteFeedbackGains(const Eigen::MatrixXd& Qaa, 
                                 Eigen::MatrixXd& Kaq, 
                                 Eigen::MatrixXd& Kav, 
                                 Eigen::MatrixXd& Kfq, 
                                 Eigen::MatrixXd& Kfv,
                                 Eigen::MatrixXd& Kmuq, 
                                 Eigen::MatrixXd& Kmuv);

  void computeFeedforwardTerms(const Eigen::MatrixXd& Qff, 
                               const Eigen::MatrixXd& Qaf);

private:
  bool has_floating_base_;
  int dimv_, max_dimf_, dim_passive_;
  Eigen::MatrixXd Qff_inv_, Saa_, Saa_inv_, Saf_, D_hat_, D_hat_inv_,
                  L_U_, L_L_, Kaq_, Kav_, Kfq_, Kfv_, Kmuq_, Kmuv_;
};

} // namespace idocp


#endif // IDOCP_BLOCK_MATRIX_INVERTER_HPP_