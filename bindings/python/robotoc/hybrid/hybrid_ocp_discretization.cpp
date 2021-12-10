#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/eigen.h>
#include <pybind11/numpy.h>

#include "robotoc/hybrid/discretization_method.hpp"
#include "robotoc/hybrid/hybrid_ocp_discretization.hpp"

#include <iostream>


namespace robotoc {
namespace python {

namespace py = pybind11;

PYBIND11_MODULE(hybrid_ocp_discretization, m) {
  py::enum_<DiscretizationMethod>(m, "DiscretizationMethod", py::arithmetic())
    .value("GridBased",  DiscretizationMethod::GridBased)
    .value("PhaseBased", DiscretizationMethod::PhaseBased)
    .export_values();

  py::class_<HybridOCPDiscretization>(m, "HybridOCPDiscretization")
    .def(py::init<const double, const int, const int>(), 
          py::arg("T"), py::arg("N"), py::arg("max_num_each_discrete_events"))
    .def("set_discretization_method", &HybridOCPDiscretization::setDiscretizationMethod,
          py::arg("discretization_method"))
    .def("discretize", &HybridOCPDiscretization::discretize,
          py::arg("contact_sequence"), py::arg("t"))
    .def("mesh_refinement", &HybridOCPDiscretization::meshRefinement,
          py::arg("contact_sequence"), py::arg("t"))
    .def("N", &HybridOCPDiscretization::N)
    .def("N_impulse", &HybridOCPDiscretization::N_impulse)
    .def("N_lift", &HybridOCPDiscretization::N_lift)
    .def("N_ideal", &HybridOCPDiscretization::N_ideal)
    .def("N_phase", &HybridOCPDiscretization::N_phase,
          py::arg("phase"))
    .def("num_contact_phases", &HybridOCPDiscretization::numContactPhases)
    .def("num_discrete_events", &HybridOCPDiscretization::numDiscreteEvents)
    .def("contact_phase", &HybridOCPDiscretization::contactPhase,
          py::arg("time_stage"))
    .def("contact_phase_after_impulse", &HybridOCPDiscretization::contactPhaseAfterImpulse,
          py::arg("impulse_index"))
    .def("contact_phase_after_lift", &HybridOCPDiscretization::contactPhaseAfterLift,
          py::arg("lift_index"))
    .def("impulse_index_after_time_stage", &HybridOCPDiscretization::impulseIndexAfterTimeStage,
          py::arg("time_stage"))
    .def("lift_index_after_time_stage", &HybridOCPDiscretization::liftIndexAfterTimeStage,
          py::arg("time_stage"))
    .def("time_stage_before_impulse", &HybridOCPDiscretization::timeStageBeforeImpulse,
          py::arg("impulse_index"))
    .def("time_stage_after_impulse", &HybridOCPDiscretization::timeStageAfterImpulse,
          py::arg("impulse_index"))
    .def("time_stage_before_lift", &HybridOCPDiscretization::timeStageBeforeLift,
          py::arg("lift_index"))
    .def("time_stage_after_lift", &HybridOCPDiscretization::timeStageAfterLift,
          py::arg("lift_index"))
    .def("is_time_stage_before_impulse", &HybridOCPDiscretization::isTimeStageBeforeImpulse,
          py::arg("time_stage"))
    .def("is_time_stage_after_impulse", &HybridOCPDiscretization::isTimeStageAfterImpulse,
          py::arg("time_stage"))
    .def("is_time_stage_before_lift", &HybridOCPDiscretization::isTimeStageBeforeLift,
          py::arg("time_stage"))
    .def("is_time_stage_after_lift", &HybridOCPDiscretization::isTimeStageAfterLift,
          py::arg("time_stage"))
    .def("t", &HybridOCPDiscretization::t,
          py::arg("time_stage"))
    .def("t_impulse", &HybridOCPDiscretization::t_impulse,
          py::arg("impulse_index"))
    .def("t_lift", &HybridOCPDiscretization::t_lift,
          py::arg("lift_index"))
    .def("dt", &HybridOCPDiscretization::dt,
          py::arg("time_stage"))
    .def("dt_aux", &HybridOCPDiscretization::dt_aux,
          py::arg("impulse_index"))
    .def("dt_lift", &HybridOCPDiscretization::dt_lift,
          py::arg("lift_index"))
    .def("dt_max", &HybridOCPDiscretization::dt_max)
    .def("dt_ideal", &HybridOCPDiscretization::dt_ideal)
    .def("is_STO_enabled_event", &HybridOCPDiscretization::isSTOEnabledEvent,
          py::arg("event_index"))
    .def("is_STO_enabled_phase", &HybridOCPDiscretization::isSTOEnabledPhase,
          py::arg("phase"))
    .def("is_STO_enabled_next_phase", &HybridOCPDiscretization::isSTOEnabledNextPhase,
          py::arg("phase"))
    .def("is_STO_enabled_stage", &HybridOCPDiscretization::isSTOEnabledStage,
          py::arg("time_stage"))
    .def("is_STO_enabled_impulse", &HybridOCPDiscretization::isSTOEnabledImpulse,
          py::arg("impulse_index"))
    .def("is_STO_enabled_lift", &HybridOCPDiscretization::isSTOEnabledLift,
          py::arg("lift_index"))
    .def("event_index_impulse", &HybridOCPDiscretization::eventIndexImpulse,
          py::arg("impulse_index"))
    .def("event_index_lift", &HybridOCPDiscretization::eventIndexLift,
          py::arg("lift_index"))
    .def("event_type", &HybridOCPDiscretization::eventType,
          py::arg("event_index"))
    .def("discretization_method", &HybridOCPDiscretization::discretizationMethod)
    .def("max_num_each_discrete_events", &HybridOCPDiscretization::maxNumEachDiscreteEvents)
    .def("time_steps", &HybridOCPDiscretization::timeSteps)
    .def("time_points", &HybridOCPDiscretization::timePoints)
    .def("is_formulation_tractable", &HybridOCPDiscretization::isFormulationTractable)
    .def("is_switching_time_consistent", &HybridOCPDiscretization::isSwitchingTimeConsistent)
    .def("__str__", [](const HybridOCPDiscretization& self) {
        std::stringstream ss;
        ss << self;
        return ss.str();
      });
}

} // namespace python
} // namespace robotoc