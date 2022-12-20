#include <pybind11/pybind11.h>

#include "robotoc/constraints/constraint_component_base.hpp"
#include "robotoc/utils/pybind11_macros.hpp"


namespace robotoc {
namespace python {

class PyConstraintComponentBase : public ConstraintComponentBase {
public:
  // Inherit the constructors
  using ConstraintComponentBase::ConstraintComponentBase;

  KinematicsLevel kinematicsLevel() const override {
    PYBIND11_OVERRIDE_PURE(KinematicsLevel, ConstraintComponentBase, 
                           kinematicsLevel, );
  }

  void allocateExtraData(ConstraintComponentData& data) const override {
    PYBIND11_OVERRIDE_PURE(void, ConstraintComponentBase, 
                           allocateExtraData, 
                           data);
  }

  bool isFeasible(Robot& robot, const ContactStatus& contact_status, 
                  const GridInfo& grid_info, const SplitSolution& s,
                  ConstraintComponentData& data) const override {
    PYBIND11_OVERRIDE_PURE(bool, ConstraintComponentBase, 
                           isFeasible, 
                           robot, contact_status, grid_info, s, data);
  }

  void setSlack(Robot& robot, const ContactStatus& contact_status, 
                const GridInfo& grid_info, const SplitSolution& s,
                ConstraintComponentData& data) const override {
    PYBIND11_OVERRIDE_PURE(void, ConstraintComponentBase, 
                           setSlack, 
                           robot, contact_status, grid_info, s, data);
  }

  void evalConstraint(Robot& robot, const ContactStatus& contact_status, 
                      const GridInfo& grid_info, const SplitSolution& s,
                      ConstraintComponentData& data) const override {
    PYBIND11_OVERRIDE_PURE(void, ConstraintComponentBase, 
                           evalConstraint, 
                           robot, contact_status, grid_info, s, data);
  }

  void evalDerivatives(Robot& robot, const ContactStatus& contact_status, 
                       const GridInfo& grid_info, const SplitSolution& s,
                       ConstraintComponentData& data,
                       SplitKKTResidual& kkt_residual) const override {
    PYBIND11_OVERRIDE_PURE(void, ConstraintComponentBase, 
                           evalDerivatives, 
                           robot, contact_status, grid_info, s, data, kkt_residual);
  }

  void condenseSlackAndDual(const ContactStatus& contact_status,
                            const GridInfo& grid_info,
                            ConstraintComponentData& data, 
                            SplitKKTMatrix& kkt_matrix,
                            SplitKKTResidual& kkt_residual) const override {
    PYBIND11_OVERRIDE_PURE(void, ConstraintComponentBase, 
                           condenseSlackAndDual, 
                           contact_status, grid_info, data, kkt_matrix, kkt_residual);
  }

  void expandSlackAndDual(const ContactStatus& contact_status, 
                          const GridInfo& grid_info,
                          const SplitDirection& d,
                          ConstraintComponentData& data) const override {
    PYBIND11_OVERRIDE_PURE(void, ConstraintComponentBase, 
                           expandSlackAndDual, 
                           contact_status, grid_info, data, d);
  }

  int dimc() const override {
    PYBIND11_OVERRIDE_PURE(int, ConstraintComponentBase, 
                           dimc, );
  }
};


namespace py = pybind11;

PYBIND11_MODULE(constraint_component_base, m) {
  py::enum_<KinematicsLevel>(m, "KinematicsLevel", py::arithmetic())
    .value("PositionLevel",  KinematicsLevel::PositionLevel)
    .value("VelocityLevel", KinematicsLevel::VelocityLevel)
    .value("AccelerationLevel", KinematicsLevel::AccelerationLevel)
    .export_values();

  py::class_<ConstraintComponentBase, 
             PyConstraintComponentBase,
             std::shared_ptr<ConstraintComponentBase>>(m, "ConstraintComponentBase")
    .def(py::init<const double, const double>(),
          py::arg("barrier_param")=1.0e-03, py::arg("fraction_to_boundary_rule")=0.995)
    .def("kinematicsLevel", &ConstraintComponentBase::kinematicsLevel)
    .def("allocateExtraData", &ConstraintComponentBase::allocateExtraData,
          py::arg("data"))
    .def("isFeasible", &ConstraintComponentBase::isFeasible,
          py::arg("robot"), py::arg("contact_status"), py::arg("grid_info"), py::arg("s"), py::arg("data"))
    .def("setSlack", &ConstraintComponentBase::setSlack,
          py::arg("robot"), py::arg("contact_status"), py::arg("grid_info"), py::arg("s"), py::arg("data"))
    .def("evalConstraint", &ConstraintComponentBase::evalConstraint,
          py::arg("robot"), py::arg("contact_status"), py::arg("grid_info"), py::arg("s"), py::arg("data"))
    .def("evalDerivatives", &ConstraintComponentBase::evalDerivatives,
          py::arg("robot"), py::arg("contact_status"), py::arg("grid_info"), py::arg("s"), py::arg("data"),
          py::arg("kkt_residual"))
    .def("condenseSlackAndDual", &ConstraintComponentBase::condenseSlackAndDual,
          py::arg("contact_status"), py::arg("grid_info"), py::arg("data"), py::arg("kkt_matrix"), py::arg("kkt_residual"))
    .def("expandSlackAndDual", &ConstraintComponentBase::expandSlackAndDual,
          py::arg("contact_status"), py::arg("grid_info"), py::arg("d"), py::arg("data"))
    .def("dimc", &ConstraintComponentBase::dimc);

}

} // namespace python
} // namespace robotoc