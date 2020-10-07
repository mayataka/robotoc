#ifndef IDOCP_MPC_HPP_
#define IDOCP_MPC_HPP_ 

#include <memory>

#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/cost/cost_function.hpp"
#include "idocp/constraints/constraints.hpp"


namespace idocp {

template <typename OCPType>
class MPC {
public:

  ///
  /// @brief Construct a MPC.
  /// @param[in] robot Robot model. Must be initialized by URDF or XML.
  /// @param[in] cost Shared ptr to the cost function.
  /// @param[in] constraints Shared ptr to the constraints.
  /// @param[in] T Length of the horizon. Must be positive.
  /// @param[in] N Number of discretization of the horizon. Must be more than 1. 
  /// @param[in] num_proc Number of the threads in solving the OCP. Must be 
  /// positive. Default is 1.
  ///
  MPC(const Robot& robot, const std::shared_ptr<CostFunction>& cost,
      const std::shared_ptr<Constraints>& constraints, const double T, 
      const int N, const int num_proc=1);

  ///
  /// @brief Default constructor. Does not construct any datas. 
  ///
  MPC();

  ///
  /// @brief Destructor. 
  ///
  ~MPC();

  ///
  /// @brief Use default copy constructor. 
  ///
  MPC(const MPC&) = default;

  ///
  /// @brief Use default copy assign operator. 
  ///
  MPC& operator=(const MPC&) = default;

  ///
  /// @brief Use default move constructor. 
  ///
  MPC(MPC&&) noexcept = default;

  ///
  /// @brief Use default move assign operator. 
  ///
  MPC& operator=(MPC&&) noexcept = default;

  ///
  /// @brief Initializes the solution. First, sets the configuration and 
  /// velocity over the horizon uniformly. If the OCP class is ParNMPC, the 
  /// auxiliary matrix is also initialized by the hessian of the terminal cost.
  /// Then the OCP with line search is solved at the predefine times.
  /// @param[in] t Current time. 
  /// @param[in] q Configuration. Size must be Robot::dimq().
  /// @param[in] v Velocity. Size must be Robot::dimv().
  /// @param[in] max_itr Iteration times of the OCP to initialize the solution. 
  /// default is 0, i.e., does not solve OCP in default setting.
  ///
  void initializeSolution(const double t, const Eigen::VectorXd& q, 
                          const Eigen::VectorXd& v, const int max_itr=0);

  ///
  /// @brief Activate a contact over specified time steps 
  /// (from time_stage_begin until time_stage_end). 
  /// @param[in] contact_index Index of a contact of interedted. 
  /// @param[in] time_stage_begin Beginning of time stages to activate. 
  /// @param[in] time_stage_end End of time stages to activate. 
  ///
  void activateContact(const int contact_index, const int time_stage_begin, 
                       const int time_stage_end);

  ///
  /// @brief Deactivate a contact over specified time steps 
  /// (from time_stage_begin until time_stage_end). 
  /// @param[in] contact_index Index of a contact of interedted. 
  /// @param[in] time_stage_begin Beginning of time stages to deactivate. 
  /// @param[in] time_stage_end End of time stages to deactivate. 
  ///
  void deactivateContact(const int contact_index, const int time_stage_begin, 
                         const int time_stage_end);

  ///
  /// @brief Activate contacts over specified time steps 
  /// (from time_stage_begin until time_stage_end). 
  /// @param[in] contact_indices Indices of contacts of interedted. 
  /// @param[in] time_stage_begin Beginning of time stages to activate. 
  /// @param[in] time_stage_end End of time stages to activate. 
  ///
  void activateContacts(const std::vector<int>& contact_indices, 
                        const int time_stage_begin, const int time_stage_end);

  ///
  /// @brief Deactivate contacts over specified time steps 
  /// (from time_stage_begin until time_stage_end). 
  /// @param[in] contact_indices Indices of contacts of interedted. 
  /// @param[in] time_stage_begin Beginning of time stages to deactivate. 
  /// @param[in] time_stage_end End of time stages to deactivate. 
  ///
  void deactivateContacts(const std::vector<int>& contact_indices, 
                          const int time_stage_begin, const int time_stage_end);

  ///
  /// @brief Sets the contact points over the horizon. 
  /// @param[in] contact_points Contact points over the horizon.
  ///
  void setContactPoint(const std::vector<Eigen::Vector3d>& contact_points);


  ///
  /// @brief Sets the contact points over the horizon by the configuration. 
  /// @param[in] q configuration. Size must be Robot::dimq().
  ///
  void setContactPointByKinematics(const Eigen::VectorXd& q);

  ///
  /// @brief Updates solution by computing the primal-dual Newon direction.
  /// @param[in] t Current time. 
  /// @param[in] q Initial configuration. Size must be Robot::dimq().
  /// @param[in] v Initial velocity. Size must be Robot::dimv().
  ///
  void updateSolution(const double t, const Eigen::VectorXd& q, 
                      const Eigen::VectorXd& v);

  ///
  /// @brief Get the contorl input torques of the initial stage.
  /// @param[out] u The control input torques. Size must be Robot::dimv().
  ///
  void getControlInput(Eigen::VectorXd& u);

  ///
  /// @brief Gets the state-feedback gain for the control input torques at the 
  /// initial stage.
  /// @param[out] Kq Gain with respec to the configuration. Size must be 
  /// Robot::dimv() x Robot::dimv().
  /// @param[out] Kv Gain with respec to the velocity. Size must be
  /// Robot::dimv() x Robot::dimv().
  ///
  void getStateFeedbackGain(Eigen::MatrixXd& Kq, Eigen::MatrixXd& Kv);

  ///
  /// @brief Returns the squared KKT error norm by using previously computed 
  /// results computed by updateSolution(). The result is not exactly the 
  /// same as the squared KKT error norm of the original OCP. The result is the
  /// squared norm of the condensed residual. However, this variables is 
  /// sufficiently close to the original KKT error norm.
  /// @param[in] t Current time. 
  /// @return The squared norm of the condensed KKT residual.
  ///
  double KKTError(const double t);

  ///
  /// @brief Computes and returns the squared KKT error norm of the OCP. 
  /// @param[in] t Current time. 
  /// @param[in] q Initial configuration. Size must be Robot::dimq().
  /// @param[in] v Initial velocity. Size must be Robot::dimv().
  /// @return The squared norm of the kKT residual.
  ///
  double computeKKTError(const double t, const Eigen::VectorXd& q, 
                         const Eigen::VectorXd& v);

private:
  OCPType ocp_;

};

} // namespace idocp 

#include "idocp/ocp/mpc.hxx"

#endif // IDOCP_MPC_HPP_