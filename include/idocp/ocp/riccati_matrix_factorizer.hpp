#ifndef IDOCP_RICCATI_MATRIX_FACTORIZER_HPP_
#define IDOCP_RICCATI_MATRIX_FACTORIZER_HPP_

#include <vector>
#include <utility>

#include "Eigen/Core"

#include "robot/robot.hpp"


namespace idocp {

class RiccatiMatrixFactorizer {
public:
  // Constructor.
  // Argments:
  //    robot: The robot model that has been already initialized.
  RiccatiMatrixFactorizer(const Robot& robot);

  // Default constructor.
  RiccatiMatrixFactorizer();

  // Destructor.
  ~RiccatiMatrixFactorizer();
 
  // Use default copy constructor.
  RiccatiMatrixFactorizer(const RiccatiMatrixFactorizer& other) = default;

  // Use default copy operator.
  RiccatiMatrixFactorizer& operator=(const RiccatiMatrixFactorizer& other) 
      = default;

  void setIntegrationSensitivities(const Robot& robot, const double dtau,
                                   const Eigen::VectorXd& q,
                                   const Eigen::VectorXd& v);

  void factorize(const double dtau, const Eigen::MatrixXd& Pqq_next, 
                 const Eigen::MatrixXd& Pqv_next, 
                 const Eigen::MatrixXd& Pvq_next, 
                 const Eigen::MatrixXd& Pvv_next, Eigen::MatrixXd& Qqq, 
                 Eigen::MatrixXd& Qqv, Eigen::MatrixXd& Qvq, 
                 Eigen::MatrixXd& Qvv);

  void factorize(const double dtau, const Eigen::MatrixXd& Pqv_next, 
                 const Eigen::MatrixXd& Pvv_next, Eigen::MatrixXd& Qqa, 
                 Eigen::MatrixXd& Qva);

  void factorize(const double dtau, const Eigen::MatrixXd& Pvv_next, 
                 Eigen::MatrixXd& Qaa);

private:
  bool has_floating_base_;
  int dimv_;
  Eigen::MatrixXd dintegrate_dq_, dintegrate_dv_;

};

} // namespace idocp


#endif // IDOCP_RICCATI_MATRIX_FACTORIZER_HPP_