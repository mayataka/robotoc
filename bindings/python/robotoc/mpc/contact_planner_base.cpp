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

  bool plan(const double t, const Eigen::VectorXd& q, const Eigen::VectorXd& v, 
            const ContactStatus& contact_status, 
            const int planning_steps) override {
    PYBIND11_OVERRIDE_PURE(bool, ContactPlannerBase, 
                           plan, 
                           t, q, v, contact_status, planning_steps);
  }

  const aligned_vector<SE3>& contactPlacements(const int step) const override {
    PYBIND11_OVERRIDE_PURE(const aligned_vector<SE3>&, ContactPlannerBase, 
                           contactPlacements, step);
  }

  const std::vector<Eigen::Vector3d>& contactPositions(const int step) const override {
    PYBIND11_OVERRIDE_PURE(const std::vector<Eigen::Vector3d>&, ContactPlannerBase, 
                           contactPositions, step);
  }

  const std::vector<Eigen::Matrix3d>& contactSurfaces(const int step) const override {
    PYBIND11_OVERRIDE_PURE(const std::vector<Eigen::Matrix3d>&, ContactPlannerBase, 
                           contactSurfaces, step);
  }

  const Eigen::Vector3d& CoM(const int step) const override {
    PYBIND11_OVERRIDE_PURE(const Eigen::Vector3d&, ContactPlannerBase, 
                           CoM, step);
  }

  const Eigen::Matrix3d& R(const int step) const override {
    PYBIND11_OVERRIDE_PURE(const Eigen::Matrix3d&, ContactPlannerBase, 
                           R, step);
  }

  int size() const override {
    PYBIND11_OVERRIDE_PURE(int, ContactPlannerBase, 
                           size, );
  }
};


namespace py = pybind11;

PYBIND11_MODULE(contact_planner_base, m) {
  py::class_<ContactPlannerBase, 
             PyContactPlannerBase, 
             std::shared_ptr<ContactPlannerBase>>(m, "ContactPlannerBase")
    .def(py::init<>()) 
    .def("plan", &ContactPlannerBase::plan,
          py::arg("t"), py::arg("q"), py::arg("v"), py::arg("contact_status"), 
          py::arg("planning_steps"))
    .def("contactPlacements", &ContactPlannerBase::contactPlacements,
          py::arg("step"))
    .def("contactPositions", &ContactPlannerBase::contactPositions,
          py::arg("step"))
    .def("contactSurfaces", &ContactPlannerBase::contactSurfaces,
          py::arg("step"))
    .def("com", &ContactPlannerBase::CoM,
          py::arg("step"))
    .def("R", &ContactPlannerBase::R,
          py::arg("step"))
    .def("size", &ContactPlannerBase::size)
    .def("__str__", [](const ContactPlannerBase& self) {
        std::stringstream ss;
        self.disp(ss);
        return ss.str();
      });
}

} // namespace python
} // namespace robotoc