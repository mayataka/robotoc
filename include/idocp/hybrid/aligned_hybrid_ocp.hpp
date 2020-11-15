#ifndef IDOCP_ALIGNED_HYBRID_OCP_HPP_
#define IDOCP_ALIGNED_HYBRID_OCP_HPP_

#include "idocp/ocp/split_ocp.hpp"
#include "idocp/impulse/split_impulse_ocp.hpp"
#include "idocp/ocp/terminal_ocp.hpp"
#include "idocp/hybrid/contact_sequence.hpp"
#include "idocp/hybrid/hybrid_container.hpp"

namespace idocp {

class AlignedHybridOCP {
public:
  HybridOCP();
  ~(HybridOCP);

  void linearizeOCP();

  void update(const ContactSequence& contact_sequence) {

  }

  int N() const;

  double dtau(const int i) const;

  double t(const int i) const;

private:
  std::vector<double> dtau_;
  std::vector<double> dtau_lift_;

};
  
} // namespace idocp

#endif // IDOCP_ALIGNED_HYBRID_OCP_HPP_ 