#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/eigen.h>
#include <pybind11/numpy.h>

#include "robotoc/mpc/contact_planner_base.hpp"


namespace robotoc {
namespace python {

class PyContactPlannerBase : public ContactPlannerBase {
public:
  // Inherit the constructors
  using ContactPlannerBase::ContactPlannerBase;

  void init(const Eigen::VectorXd& q) override {
    PYBIND11_OVERRIDE_PURE(void, ContactPlannerBase, 
                           init, q);
  }

  bool plan(const Eigen::VectorXd& q, const Eigen::VectorXd& v, 
            const ContactStatus& contact_status, 
            const int planning_steps) override {
    PYBIND11_OVERRIDE_PURE(bool, ContactPlannerBase, 
                           plan, 
                           q, v, contact_status, planning_steps);
  }

  const aligned_vector<SE3>& contactPlacement(const int step) const override {
    PYBIND11_OVERRIDE_PURE(const aligned_vector<SE3>&, ContactPlannerBase, 
                           contactPlacement, step);
  }

  const aligned_vector<aligned_vector<SE3>>& contactPlacement() const override {
    PYBIND11_OVERRIDE_PURE(const aligned_vector<aligned_vector<SE3>>&, ContactPlannerBase, 
                           contactPlacement, );
  }

  const std::vector<Eigen::Vector3d>& contactPosition(const int step) const override {
    PYBIND11_OVERRIDE_PURE(const std::vector<Eigen::Vector3d>&, ContactPlannerBase, 
                           contactPosition, step);
  }

  const std::vector<std::vector<Eigen::Vector3d>>& contactPosition() const override {
    PYBIND11_OVERRIDE_PURE(const std::vector<std::vector<Eigen::Vector3d>>&, ContactPlannerBase, 
                           contactPosition, );
  }

  const Eigen::Vector3d& com(const int step) const override {
    PYBIND11_OVERRIDE_PURE(const Eigen::Vector3d&, ContactPlannerBase, 
                           com, step);
  }

  const std::vector<Eigen::Vector3d>& com() const override {
    PYBIND11_OVERRIDE_PURE(const std::vector<Eigen::Vector3d>&, ContactPlannerBase, 
                           com, );
  }

  const Eigen::Matrix3d& R(const int step) const override {
    PYBIND11_OVERRIDE_PURE(const Eigen::Matrix3d&, ContactPlannerBase, 
                           R, step);
  }

  const std::vector<Eigen::Matrix3d>& R() const override {
    PYBIND11_OVERRIDE_PURE(const std::vector<Eigen::Matrix3d>&, ContactPlannerBase, 
                           R, );
  }
};


namespace py = pybind11;

PYBIND11_MODULE(contact_planner_base, m) {
  py::class_<ContactPlannerBase, 
             PyContactPlannerBase, 
             std::shared_ptr<ContactPlannerBase>>(m, "ContactPlannerBase")
    .def(py::init<>());
}

} // namespace python
} // namespace robotoc