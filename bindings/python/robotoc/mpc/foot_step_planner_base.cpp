#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/eigen.h>
#include <pybind11/numpy.h>

#include "robotoc/mpc/foot_step_planner_base.hpp"


namespace robotoc {
namespace python {

class PyFootStepPlannerBase : public FootStepPlannerBase {
public:
  // Inherit the constructors
  using FootStepPlannerBase::FootStepPlannerBase;

  void init(const Eigen::VectorXd& q) override {
    PYBIND11_OVERRIDE_PURE(void, FootStepPlannerBase, 
                           init, q);
  }

  bool plan(const Eigen::VectorXd& q, const ContactStatus& contact_status, 
            const int planning_steps) override {
    PYBIND11_OVERRIDE_PURE(bool, FootStepPlannerBase, 
                           plan, 
                           q, contact_status, planning_steps);
  }

  const aligned_vector<SE3>& contactPlacement(const int step) const override {
    PYBIND11_OVERRIDE_PURE(const aligned_vector<SE3>&, FootStepPlannerBase, 
                           contactPlacement, step);
  }

  const std::vector<Eigen::Vector3d>& contactPosition(const int step) const override {
    PYBIND11_OVERRIDE_PURE(const std::vector<Eigen::Vector3d>&, FootStepPlannerBase, 
                           contactPosition, step);
  }

  const Eigen::Vector3d& com(const int step) const override {
    PYBIND11_OVERRIDE_PURE(const Eigen::Vector3d&, FootStepPlannerBase, 
                           com, step);
  }
};


namespace py = pybind11;

PYBIND11_MODULE(foot_step_planner_base, m) {
  py::class_<FootStepPlannerBase, 
             PyFootStepPlannerBase, 
             std::shared_ptr<FootStepPlannerBase>>(m, "FootStepPlannerBase")
    .def(py::init<>());
}

} // namespace python
} // namespace robotoc