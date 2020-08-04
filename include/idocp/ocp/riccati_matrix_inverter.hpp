#ifndef IDOCP_BLOCK_MATRIX_INVERTER_HPP_
#define IDOCP_BLOCK_MATRIX_INVERTER_HPP_

#include <vector>
#include <utility>

#include "Eigen/Core"

#include "idocp/robot/robot.hpp"


namespace idocp {

class RiccatiMatrixInverter {
public:
  // Constructor.
  // Argments:
  //    robot: The robot model that has been already initialized.
  RiccatiMatrixInverter(const Robot& robot);

  // Default constructor.
  RiccatiMatrixInverter();

  // Destructor.
  ~RiccatiMatrixInverter();
 
  // Use default copy constructor.
  RiccatiMatrixInverter(const RiccatiMatrixInverter& other) = default;

  // Use default copy operator.
  RiccatiMatrixInverter& operator=(const RiccatiMatrixInverter& other) 
      = default;

  void setContactStatus(const Robot& robot);

  void precompute(const Eigen::MatrixXd& Qaf, const Eigen::MatrixXd& Qff);

  // Fixed base without contacts
  void invert(const Eigen::MatrixXd& Qqa, const Eigen::MatrixXd& Qva, 
              const Eigen::MatrixXd& Qaa, const Eigen::VectorXd& la,
              Eigen::MatrixXd& Kaq, Eigen::MatrixXd& Kav, Eigen::VectorXd& ka);

  // Fixed base with contacts
  void invert(const Eigen::MatrixXd& Qqa, const Eigen::MatrixXd& Qva, 
              const Eigen::MatrixXd& Qaa, const Eigen::MatrixXd& Qqf, 
              const Eigen::MatrixXd& Qvf, const Eigen::MatrixXd& Cq, 
              const Eigen::MatrixXd& Cv, const Eigen::MatrixXd& Ca, 
              const Eigen::VectorXd& la, const Eigen::VectorXd& lf, 
              const Eigen::VectorXd& C_res, Eigen::MatrixXd& Kaq, 
              Eigen::MatrixXd& Kav, Eigen::MatrixXd& Kfq, Eigen::MatrixXd& Kfv, 
              Eigen::MatrixXd& Kmuq, Eigen::MatrixXd& Kmuv, Eigen::VectorXd& ka, 
              Eigen::VectorXd& kf, Eigen::VectorXd& kmu);

  // Floating base without contacts
  void invert(const Eigen::MatrixXd& Qqa, const Eigen::MatrixXd& Qva, 
              const Eigen::MatrixXd& Qaa, 
              const Eigen::MatrixXd& Cq, 
              const Eigen::MatrixXd& Cv, 
              const Eigen::MatrixXd& Ca, 
              const Eigen::VectorXd& la, const Eigen::VectorXd& C_res, 
              Eigen::MatrixXd& Kaq, Eigen::MatrixXd& Kav, Eigen::MatrixXd& Kmuq, 
              Eigen::MatrixXd& Kmuv, Eigen::VectorXd& ka, Eigen::VectorXd& kmu);

  // Floating base with contacts
  void invert(const Eigen::MatrixXd& Qqa, const Eigen::MatrixXd& Qva, 
              const Eigen::MatrixXd& Qaa, const Eigen::MatrixXd& Qqf, 
              const Eigen::MatrixXd& Qvf, const Eigen::MatrixXd& Cq, 
              const Eigen::MatrixXd& Cv, const Eigen::MatrixXd& Ca, 
              const Eigen::MatrixXd& Cf, const Eigen::VectorXd& la, 
              const Eigen::VectorXd& lf, const Eigen::VectorXd& C_res, 
              Eigen::MatrixXd& Kaq, Eigen::MatrixXd& Kav, Eigen::MatrixXd& Kfq, 
              Eigen::MatrixXd& Kfv, Eigen::MatrixXd& Kmuq, 
              Eigen::MatrixXd& Kmuv, Eigen::VectorXd& ka, Eigen::VectorXd& kf, 
              Eigen::VectorXd& kmu);

private:
  bool has_floating_base_;
  int dimv_, max_dimf_, dimf_, dim_passive_;
  Eigen::MatrixXd Sff_inv_, Saa_, Saa_inv_, Saf_, Sac_, Sfc_, Scc_, Scc_inv_;
};

} // namespace idocp


#endif // IDOCP_BLOCK_MATRIX_INVERTER_HPP_