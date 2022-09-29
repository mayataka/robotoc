#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/eigen.h>
#include <pybind11/numpy.h>

#include "robotoc/ocp/discretization_method.hpp"
#include "robotoc/ocp/time_discretization.hpp"

#include <iostream>


namespace robotoc {
namespace python {

namespace py = pybind11;

PYBIND11_MODULE(time_discretization, m) {
  py::enum_<DiscretizationMethod>(m, "DiscretizationMethod", py::arithmetic())
    .value("GridBased",  DiscretizationMethod::GridBased)
    .value("PhaseBased", DiscretizationMethod::PhaseBased)
    .export_values();

  py::class_<TimeDiscretization>(m, "TimeDiscretization")
    .def(py::init<const double, const int, const int>(), 
          py::arg("T"), py::arg("N"), py::arg("reserved_num_discrete_events")=0)
    .def("set_discretization_method", &TimeDiscretization::setDiscretizationMethod,
          py::arg("discretization_method"))
    .def("discretize", &TimeDiscretization::discretize,
          py::arg("contact_sequence"), py::arg("t"))
    .def("mesh_refinement", &TimeDiscretization::meshRefinement,
          py::arg("contact_sequence"), py::arg("t"))
    .def("N", &TimeDiscretization::N)
    .def("N_impulse", &TimeDiscretization::N_impulse)
    .def("N_lift", &TimeDiscretization::N_lift)
    .def("N_ideal", &TimeDiscretization::N_ideal)
    .def("N_phase", &TimeDiscretization::N_phase,
          py::arg("phase"))
    .def("num_contact_phases", &TimeDiscretization::numContactPhases)
    .def("num_discrete_events", &TimeDiscretization::numDiscreteEvents)
    .def("contact_phase", &TimeDiscretization::contactPhase,
          py::arg("time_stage"))
    .def("contact_phase_after_impulse", &TimeDiscretization::contactPhaseAfterImpulse,
          py::arg("impulse_index"))
    .def("contact_phase_after_lift", &TimeDiscretization::contactPhaseAfterLift,
          py::arg("lift_index"))
    .def("impulse_index_after_time_stage", &TimeDiscretization::impulseIndexAfterTimeStage,
          py::arg("time_stage"))
    .def("lift_index_after_time_stage", &TimeDiscretization::liftIndexAfterTimeStage,
          py::arg("time_stage"))
    .def("time_stage_before_impulse", &TimeDiscretization::timeStageBeforeImpulse,
          py::arg("impulse_index"))
    .def("time_stage_after_impulse", &TimeDiscretization::timeStageAfterImpulse,
          py::arg("impulse_index"))
    .def("time_stage_before_lift", &TimeDiscretization::timeStageBeforeLift,
          py::arg("lift_index"))
    .def("time_stage_after_lift", &TimeDiscretization::timeStageAfterLift,
          py::arg("lift_index"))
    .def("is_time_stage_before_impulse", &TimeDiscretization::isTimeStageBeforeImpulse,
          py::arg("time_stage"))
    .def("is_time_stage_after_impulse", &TimeDiscretization::isTimeStageAfterImpulse,
          py::arg("time_stage"))
    .def("is_time_stage_before_lift", &TimeDiscretization::isTimeStageBeforeLift,
          py::arg("time_stage"))
    .def("is_time_stage_after_lift", &TimeDiscretization::isTimeStageAfterLift,
          py::arg("time_stage"))
    .def("t0", &TimeDiscretization::t0)
    .def("tf", &TimeDiscretization::tf)
    .def("impulseTime", &TimeDiscretization::impulseTime,
          py::arg("impulse_index"))
    .def("liftTime", &TimeDiscretization::liftTime,
          py::arg("lift_index"))
    .def("dt_max", &TimeDiscretization::dt_max)
    .def("dt_ideal", &TimeDiscretization::dt_ideal)
    .def("grid_info", &TimeDiscretization::gridInfo,
          py::arg("time_stage"))
    .def("grid_info_impulse", &TimeDiscretization::gridInfoImpulse,
          py::arg("impulse_index"))
    .def("grid_info_aux", &TimeDiscretization::gridInfoAux,
          py::arg("aux_index"))
    .def("grid_info_lift", &TimeDiscretization::gridInfoLift,
          py::arg("lift_index"))
    .def("is_STO_enabled_event", &TimeDiscretization::isSTOEnabledEvent,
          py::arg("event_index"))
    .def("is_STO_enabled_phase", &TimeDiscretization::isSTOEnabledPhase,
          py::arg("phase"))
    .def("is_STO_enabled_next_phase", &TimeDiscretization::isSTOEnabledNextPhase,
          py::arg("phase"))
    .def("is_STO_enabled_impulse", &TimeDiscretization::isSTOEnabledImpulse,
          py::arg("impulse_index"))
    .def("is_STO_enabled_lift", &TimeDiscretization::isSTOEnabledLift,
          py::arg("lift_index"))
    .def("event_index_impulse", &TimeDiscretization::eventIndexImpulse,
          py::arg("impulse_index"))
    .def("event_index_lift", &TimeDiscretization::eventIndexLift,
          py::arg("lift_index"))
    .def("event_type", &TimeDiscretization::eventType,
          py::arg("event_index"))
    .def("discretization_method", &TimeDiscretization::discretizationMethod)
    .def("reserve", &TimeDiscretization::reserve,
          py::arg("reserved_num_discrete_events"))
    .def("reserved_num_discrete_events", &TimeDiscretization::reservedNumDiscreteEvents)
    .def("time_steps", &TimeDiscretization::timeSteps)
    .def("time_points", &TimeDiscretization::timePoints)
    .def("is_formulation_tractable", &TimeDiscretization::isFormulationTractable)
    .def("is_switching_time_consistent", &TimeDiscretization::isSwitchingTimeConsistent)
    .def("discretizeGrid", &TimeDiscretization::discretizeGrid,
          py::arg("contact_sequence"), py::arg("t")) 
    .def("getGrid", &TimeDiscretization::getGrid) 
    .def("__str__", [](const TimeDiscretization& self) {
        std::stringstream ss;
        ss << self;
        return ss.str();
      });
}

} // namespace python
} // namespace robotoc