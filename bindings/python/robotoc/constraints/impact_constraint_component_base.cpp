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
                  const GridInfo& grid_info, const SplitSolution& s,
                  ConstraintComponentData& data) const override {
    PYBIND11_OVERRIDE_PURE(bool, ImpactConstraintComponentBase, 
                           isFeasible, 
                           robot, impact_status, grid_info, s, data);
  }

  void setSlack(Robot& robot, const ImpactStatus& impact_status, 
                const GridInfo& grid_info, const SplitSolution& s,
                ConstraintComponentData& data) const override {
    PYBIND11_OVERRIDE_PURE(void, ImpactConstraintComponentBase, 
                           setSlack, 
                           robot, impact_status, grid_info, s, data);
  }

  void evalConstraint(Robot& robot, const ImpactStatus& impact_status, 
                      const GridInfo& grid_info, const SplitSolution& s,
                      ConstraintComponentData& data) const override {
    PYBIND11_OVERRIDE_PURE(void, ImpactConstraintComponentBase, 
                           evalConstraint, 
                           robot, impact_status, grid_info, s, data);
  }

  void evalDerivatives(Robot& robot, const ImpactStatus& impact_status, 
                      const GridInfo& grid_info, const SplitSolution& s,
                      ConstraintComponentData& data,
                      SplitKKTResidual& kkt_residual) const override {
    PYBIND11_OVERRIDE_PURE(void, ImpactConstraintComponentBase, 
                           evalDerivatives, 
                           robot, impact_status, grid_info, s, data, kkt_residual);
  }

  void condenseSlackAndDual(const ImpactStatus& impact_status, 
                            const GridInfo& grid_info, 
                            ConstraintComponentData& data,
                            SplitKKTMatrix& kkt_matrix,
                            SplitKKTResidual& kkt_residual) const override {
    PYBIND11_OVERRIDE_PURE(void, ImpactConstraintComponentBase, 
                           condenseSlackAndDual, 
                           impact_status, grid_info, data, kkt_matrix, kkt_residual);
  }

  void expandSlackAndDual(const ImpactStatus& impact_status, 
                          const GridInfo& grid_info, 
                          const SplitDirection& d,
                          ConstraintComponentData& data) const override {
    PYBIND11_OVERRIDE_PURE(void, ImpactConstraintComponentBase, 
                           expandSlackAndDual, 
                           impact_status, grid_info, d, data);
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
          py::arg("robot"), py::arg("contact_status"), py::arg("grid_info"), py::arg("s"), py::arg("data"))
    .def("setSlack", &ImpactConstraintComponentBase::setSlack,
          py::arg("robot"), py::arg("contact_status"), py::arg("grid_info"), py::arg("s"), py::arg("data"))
    .def("evalConstraint", &ImpactConstraintComponentBase::evalConstraint,
          py::arg("robot"), py::arg("contact_status"), py::arg("grid_info"), py::arg("s"), py::arg("data"))
    .def("evalDerivatives", &ImpactConstraintComponentBase::evalDerivatives,
          py::arg("robot"), py::arg("contact_status"), py::arg("grid_info"), py::arg("s"), py::arg("data"),
          py::arg("kkt_residual"))
    .def("condenseSlackAndDual", &ImpactConstraintComponentBase::condenseSlackAndDual,
          py::arg("contact_status"), py::arg("grid_info"), py::arg("data"), py::arg("kkt_matrix"), py::arg("kkt_residual"))
    .def("expandSlackAndDual", &ImpactConstraintComponentBase::expandSlackAndDual,
          py::arg("contact_status"), py::arg("grid_info"), py::arg("d"), py::arg("data"))
    .def("dimc", &ImpactConstraintComponentBase::dimc);
}

} // namespace python
} // namespace robotoc