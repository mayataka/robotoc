#include <pybind11/pybind11.h>

#include "robotoc/constraints/impact_constraint_component_base.hpp"
#include "robotoc/utils/pybind11_macros.hpp"


namespace robotoc {
namespace python {

class PyImpactConstraintComponentBase : public ImpactConstraintComponentBase {
public:
  // Inherit the constructors
  using ImpactConstraintComponentBase::ImpactConstraintComponentBase;

  KinematicsLevel kinematicsLevel() const override {
    PYBIND11_OVERRIDE_PURE(KinematicsLevel, ImpactConstraintComponentBase, 
                           kinematicsLevel, );
  }

  void allocateExtraData(ConstraintComponentData& data) const override {
    PYBIND11_OVERRIDE_PURE(void, ImpactConstraintComponentBase, 
                           allocateExtraData, data);
  }

  bool isFeasible(Robot& robot, const ImpactStatus& impact_status, 
                  ConstraintComponentData& data, 
                  const SplitSolution& s) const override {
    PYBIND11_OVERRIDE_PURE(bool, ImpactConstraintComponentBase, 
                           isFeasible, 
                           robot, impact_status, data, s);
  }

  void setSlack(Robot& robot, const ImpactStatus& impact_status, 
                ConstraintComponentData& data, 
                const SplitSolution& s) const override {
    PYBIND11_OVERRIDE_PURE(void, ImpactConstraintComponentBase, 
                           setSlack, 
                           robot, impact_status, data, s);
  }

  void evalConstraint(Robot& robot, const ImpactStatus& impact_status, 
                      ConstraintComponentData& data, 
                      const SplitSolution& s) const override {
    PYBIND11_OVERRIDE_PURE(void, ImpactConstraintComponentBase, 
                           evalConstraint, 
                           robot, impact_status, data, s);
  }

  void evalDerivatives(Robot& robot, const ImpactStatus& impact_status, 
                       ConstraintComponentData& data,
                       const SplitSolution& s,
                       SplitKKTResidual& kkt_residual) const override {
    PYBIND11_OVERRIDE_PURE(void, ImpactConstraintComponentBase, 
                           evalDerivatives, 
                           robot, impact_status, data, s, kkt_residual);
  }

  void condenseSlackAndDual(const ImpactStatus& impact_status, 
                            ConstraintComponentData& data,
                            SplitKKTMatrix& kkt_matrix,
                            SplitKKTResidual& kkt_residual) const override {
    PYBIND11_OVERRIDE_PURE(void, ImpactConstraintComponentBase, 
                           condenseSlackAndDual, 
                           impact_status, data, kkt_matrix, kkt_residual);
  }

  void expandSlackAndDual(const ImpactStatus& impact_status, 
                          ConstraintComponentData& data, 
                          const SplitDirection& d) const override {
    PYBIND11_OVERRIDE_PURE(void, ImpactConstraintComponentBase, 
                           expandSlackAndDual, 
                           impact_status, data, d);
  }

  int dimc() const override {
    PYBIND11_OVERRIDE_PURE(int, ImpactConstraintComponentBase, 
                           dimc, );
  }
};


namespace py = pybind11;

PYBIND11_MODULE(impact_constraint_component_base, m) {
  py::class_<ImpactConstraintComponentBase, 
             PyImpactConstraintComponentBase, 
             std::shared_ptr<ImpactConstraintComponentBase>>(m, "ImpactConstraintComponentBase")
    .def(py::init<const double, const double>(),
          py::arg("barrier_param")=1.0e-03, py::arg("fraction_to_boundary_rule")=0.995)
    .def("kinematicsLevel", &ImpactConstraintComponentBase::kinematicsLevel)
    .def("allocateExtraData", &ImpactConstraintComponentBase::allocateExtraData,
          py::arg("data"))
    .def("isFeasible", &ImpactConstraintComponentBase::isFeasible,
          py::arg("robot"), py::arg("contact_status"), py::arg("data"), py::arg("s"))
    .def("setSlack", &ImpactConstraintComponentBase::setSlack,
          py::arg("robot"), py::arg("contact_status"), py::arg("data"), py::arg("s"))
    .def("evalConstraint", &ImpactConstraintComponentBase::evalConstraint,
          py::arg("robot"), py::arg("contact_status"), py::arg("data"), py::arg("s"))
    .def("evalDerivatives", &ImpactConstraintComponentBase::evalDerivatives,
          py::arg("robot"), py::arg("contact_status"), py::arg("data"), py::arg("s"),
          py::arg("kkt_residual"))
    .def("condenseSlackAndDual", &ImpactConstraintComponentBase::condenseSlackAndDual,
          py::arg("contact_status"), py::arg("data"), py::arg("kkt_matrix"), 
          py::arg("kkt_residual"))
    .def("expandSlackAndDual", &ImpactConstraintComponentBase::expandSlackAndDual,
          py::arg("contact_status"), py::arg("data"), py::arg("d"))
    .def("dimc", &ImpactConstraintComponentBase::dimc);
}

} // namespace python
} // namespace robotoc