#include "idocp/utils/raisim_wrapper.hpp"

#include <assert.h>


namespace idocp {
namespace raisimwrapper {

void pino2rai(const Eigen::VectorXd& q_pinocchio, 
              const Eigen::VectorXd& v_pinocchio, Eigen::VectorXd& q_raisim, 
              Eigen::VectorXd& v_raisim) {
  assert(q_pinocchio.size() == 19);
  assert(v_pinocchio.size() == 18);
  assert(q_raisim.size() == 19);
  assert(v_raisim.size() == 18);
  q_raisim.coeffRef(0)  = q_pinocchio.coeff(0);
  q_raisim.coeffRef(1)  = q_pinocchio.coeff(1);
  q_raisim.coeffRef(2)  = q_pinocchio.coeff(2);
  q_raisim.coeffRef(3)  = q_pinocchio.coeff(6);
  q_raisim.coeffRef(4)  = q_pinocchio.coeff(3);
  q_raisim.coeffRef(5)  = q_pinocchio.coeff(4);
  q_raisim.coeffRef(6)  = q_pinocchio.coeff(5);
  q_raisim.coeffRef(7)  = q_pinocchio.coeff(7);
  q_raisim.coeffRef(8)  = q_pinocchio.coeff(8);
  q_raisim.coeffRef(9)  = q_pinocchio.coeff(9);
  q_raisim.coeffRef(10) = q_pinocchio.coeff(13);
  q_raisim.coeffRef(11) = q_pinocchio.coeff(14);
  q_raisim.coeffRef(12) = q_pinocchio.coeff(15);
  q_raisim.coeffRef(13) = q_pinocchio.coeff(10);
  q_raisim.coeffRef(14) = q_pinocchio.coeff(11);
  q_raisim.coeffRef(15) = q_pinocchio.coeff(12);
  q_raisim.coeffRef(16) = q_pinocchio.coeff(16);
  q_raisim.coeffRef(17) = q_pinocchio.coeff(17);
  q_raisim.coeffRef(18) = q_pinocchio.coeff(18);
  v_raisim.coeffRef(0)  = v_pinocchio.coeff(0);
  v_raisim.coeffRef(1)  = v_pinocchio.coeff(1);
  v_raisim.coeffRef(2)  = v_pinocchio.coeff(2);
  v_raisim.coeffRef(3)  = v_pinocchio.coeff(3);
  v_raisim.coeffRef(4)  = v_pinocchio.coeff(4);
  v_raisim.coeffRef(5)  = v_pinocchio.coeff(5);
  v_raisim.coeffRef(6)  = v_pinocchio.coeff(6);
  v_raisim.coeffRef(7)  = v_pinocchio.coeff(7);
  v_raisim.coeffRef(8)  = v_pinocchio.coeff(8);
  v_raisim.coeffRef(9)  = v_pinocchio.coeff(12);
  v_raisim.coeffRef(10) = v_pinocchio.coeff(13);
  v_raisim.coeffRef(11) = v_pinocchio.coeff(14);
  v_raisim.coeffRef(12) = v_pinocchio.coeff(9);
  v_raisim.coeffRef(13) = v_pinocchio.coeff(10);
  v_raisim.coeffRef(14) = v_pinocchio.coeff(11);
  v_raisim.coeffRef(15) = v_pinocchio.coeff(15);
  v_raisim.coeffRef(16) = v_pinocchio.coeff(16);
  v_raisim.coeffRef(17) = v_pinocchio.coeff(17);
}


void pino2rai(const Eigen::VectorXd& u_pinocchio, Eigen::VectorXd& u_raisim) {
  assert(u_pinocchio.size() == 19);
  assert(u_raisim.size() == 18);
  u_raisim.template head<6>().setZero();
  u_raisim.coeffRef(6)  = u_pinocchio.coeff(6);
  u_raisim.coeffRef(7)  = u_pinocchio.coeff(7);
  u_raisim.coeffRef(8)  = u_pinocchio.coeff(8);
  u_raisim.coeffRef(9)  = u_pinocchio.coeff(12);
  u_raisim.coeffRef(10) = u_pinocchio.coeff(13);
  u_raisim.coeffRef(11) = u_pinocchio.coeff(14);
  u_raisim.coeffRef(12) = u_pinocchio.coeff(9);
  u_raisim.coeffRef(13) = u_pinocchio.coeff(10);
  u_raisim.coeffRef(14) = u_pinocchio.coeff(11);
  u_raisim.coeffRef(15) = u_pinocchio.coeff(15);
  u_raisim.coeffRef(16) = u_pinocchio.coeff(16);
  u_raisim.coeffRef(17) = u_pinocchio.coeff(17);
}


void rai2pino(const Eigen::VectorXd& q_raisim, const Eigen::VectorXd& v_raisim, 
              Eigen::VectorXd& q_pinocchio, Eigen::VectorXd& v_pinocchio) {
  assert(q_raisim.size() == 19);
  assert(v_raisim.size() == 18);
  assert(q_pinocchio.size() == 19);
  assert(v_pinocchio.size() == 18);
  q_pinocchio.coeffRef(0)  = q_raisim.coeff(0);
  q_pinocchio.coeffRef(1)  = q_raisim.coeff(1);
  q_pinocchio.coeffRef(2)  = q_raisim.coeff(2);
  q_pinocchio.coeffRef(6)  = q_raisim.coeff(3);
  q_pinocchio.coeffRef(3)  = q_raisim.coeff(4);
  q_pinocchio.coeffRef(4)  = q_raisim.coeff(5);
  q_pinocchio.coeffRef(5)  = q_raisim.coeff(6);
  q_pinocchio.coeffRef(7)  = q_raisim.coeff(7);
  q_pinocchio.coeffRef(8)  = q_raisim.coeff(8);
  q_pinocchio.coeffRef(9)  = q_raisim.coeff(9);
  q_pinocchio.coeffRef(13) = q_raisim.coeff(10);
  q_pinocchio.coeffRef(14) = q_raisim.coeff(11);
  q_pinocchio.coeffRef(15) = q_raisim.coeff(12);
  q_pinocchio.coeffRef(10) = q_raisim.coeff(13);
  q_pinocchio.coeffRef(11) = q_raisim.coeff(14);
  q_pinocchio.coeffRef(12) = q_raisim.coeff(15);
  q_pinocchio.coeffRef(16) = q_raisim.coeff(16);
  q_pinocchio.coeffRef(17) = q_raisim.coeff(17);
  q_pinocchio.coeffRef(18) = q_raisim.coeff(18);
  v_pinocchio.coeffRef(0)  = v_raisim.coeff(0);
  v_pinocchio.coeffRef(1)  = v_raisim.coeff(1);
  v_pinocchio.coeffRef(2)  = v_raisim.coeff(2);
  v_pinocchio.coeffRef(3)  = v_raisim.coeff(3);
  v_pinocchio.coeffRef(4)  = v_raisim.coeff(4);
  v_pinocchio.coeffRef(5)  = v_raisim.coeff(5);
  v_pinocchio.coeffRef(6)  = v_raisim.coeff(6);
  v_pinocchio.coeffRef(7)  = v_raisim.coeff(7);
  v_pinocchio.coeffRef(8)  = v_raisim.coeff(8);
  v_pinocchio.coeffRef(9)  = v_raisim.coeff(12);
  v_pinocchio.coeffRef(10) = v_raisim.coeff(13);
  v_pinocchio.coeffRef(11) = v_raisim.coeff(14);
  v_pinocchio.coeffRef(12) = v_raisim.coeff(9);
  v_pinocchio.coeffRef(13) = v_raisim.coeff(10);
  v_pinocchio.coeffRef(14) = v_raisim.coeff(11);
  v_pinocchio.coeffRef(15) = v_raisim.coeff(15);
  v_pinocchio.coeffRef(16) = v_raisim.coeff(16);
  v_pinocchio.coeffRef(17) = v_raisim.coeff(17);
}

} // namespace raisimwrapper
} // namespace idocp