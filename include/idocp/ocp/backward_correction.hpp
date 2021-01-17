#ifndef IDOCP_BACKWARD_CORRECTION_HPP_
#define IDOCP_BACKWARD_CORRECTION_HPP_

#include <vector>

#include "Eigen/Core"
#include "Eigen/LU"

#include "idocp/robot/robot.hpp"
#include "idocp/ocp/split_kkt_matrix.hpp"
#include "idocp/ocp/split_kkt_residual.hpp"
#include "idocp/ocp/split_backward_correction.hpp"
#include "idocp/hybrid/hybrid_container.hpp"


namespace idocp {

///
/// @class BackwardCorrection
/// @brief Backward correction.
///
class BackwardCorrection {
public:
  ///
  /// @brief Construct factorizer.
  /// @param[in] robot Robot model. Must be initialized by URDF or XML.
  /// @param[in] N Number of discretization of the horizon. Must be more than 1. 
  /// @param[in] nthreads Number of the threads in solving the optimal control 
  /// problem. Must be positive. Default is 1.
  ///
  BackwardCorrection(const Robot& robot, const double T, const int N, 
                     const int max_num_impulse, const int nthreads);

  ///
  /// @brief Default constructor. 
  ///
  BackwardCorrection();

  ///
  /// @brief Destructor. 
  ///
  ~BackwardCorrection();
 
  ///
  /// @brief Default copy constructor. 
  ///
  BackwardCorrection(const BackwardCorrection&) = default;

  ///
  /// @brief Default copy operator. 
  ///
  BackwardCorrection& operator=(const BackwardCorrection&) = default;

  ///
  /// @brief Default move constructor. 
  ///
  BackwardCorrection(BackwardCorrection&&) noexcept = default;

  ///
  /// @brief Default move assign operator. 
  ///
  BackwardCorrection& operator=(BackwardCorrection&&) noexcept = default;

  void initAuxMat(std::vector<Robot>& robots, ParNMPC& parnmpc, 
                  const double t, const Solution& s, KKTMatrix& kkt_matrix);

  void coarseUpdate(std::vector<Robot>& robots, 
                    const ContactSequence& contact_sequence, ParNMPC& parnmpc, 
                    const double t, const Eigen::VectorXd& q, 
                    const Eigen::VectorXd& v, KKTMatrix& kkt_matrix, 
                    KKTResidual& kkt_residual,
                    const Solution& s, Direction& d);

  void backwardCorrection(std::vector<Robot>& robots, ParNMPC& parnmpc, 
                          const Solution& s, Direction& d);

  double primalStepSize() const;

  double dualStepSize() const;

private:
  int N_, max_num_impulse_, nthreads_;
  double T_, dtau_;
  BackwardCorrector corrector_;
  Solution s_new_;
  std::vector<Eigen::MatrixXd> aux_mat_, aux_mat_impulse_, aux_mat_aux_, 
                               aux_mat_lift_;
  Eigen::VectorXd primal_step_sizes_, dual_step_sizes_;

};

} // namespace idocp

#endif // IDOCP_BACKWARD_CORRECTION_HPP_ 