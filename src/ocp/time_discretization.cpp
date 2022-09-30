#include "robotoc/ocp/time_discretization.hpp"

#include <iomanip>


namespace robotoc {

void TimeDiscretization::discretizeGrid(
    const std::shared_ptr<ContactSequence>& contact_sequence, const double t) {
  const int N = N_ + contact_sequence->numLiftEvents() + 2 * contact_sequence->numImpulseEvents() + 1;
  if (grid_.size() <=N) {
    grid_.resize(N);
  }
  int next_impulse_index = 0;
  int next_lift_index = 0;
  while (next_impulse_index<contact_sequence->numImpulseEvents()) {
    if (contact_sequence->impulseTime(next_impulse_index) > t) break;
    ++next_impulse_index;
  }
  while (next_lift_index<contact_sequence->numLiftEvents()) {
    if (contact_sequence->liftTime(next_lift_index) > t) break;
    ++next_lift_index;
  }
  const double dt = T_ / N_;
  const double eps = std::sqrt(std::numeric_limits<double>::epsilon());
  const double margin = 0.5 * dt;
  int stage = 0;
  double ti = t;
  while (ti+eps < t+T_) {
    const bool has_next_impulse = (next_impulse_index < contact_sequence->numImpulseEvents());
    const bool has_next_lift = (next_lift_index < contact_sequence->numLiftEvents());
    grid_[stage].t = ti;
    grid_[stage].dt = dt;
    grid_[stage].time_stage = stage;
    grid_[stage].contact_phase = next_impulse_index + next_lift_index;
    grid_[stage].impulse_index = next_impulse_index - 1;
    grid_[stage].lift_index = next_lift_index - 1;
    grid_[stage].type = GridType::Intermediate;
    if (has_next_impulse) {
      const double next_impulse_time = contact_sequence->impulseTime(next_impulse_index);
      if (next_impulse_time <= ti+dt+eps && (next_impulse_time+margin < t+T_)) {
        grid_[stage].dt = next_impulse_time - ti;
        ++stage;
        ++next_impulse_index;
        grid_[stage].t = next_impulse_time;
        grid_[stage].dt = 0;
        grid_[stage].time_stage = stage;
        grid_[stage].contact_phase = next_impulse_index + next_lift_index;
        grid_[stage].impulse_index = next_impulse_index - 1;
        grid_[stage].lift_index = next_lift_index - 1;
        grid_[stage].type = GridType::Impulse;
        ++stage;
        grid_[stage].t = next_impulse_time;
        grid_[stage].dt = std::min(ti+dt, t+T_) - next_impulse_time;
        grid_[stage].time_stage = stage;
        grid_[stage].contact_phase = next_impulse_index + next_lift_index;
        grid_[stage].impulse_index = next_impulse_index - 1;
        grid_[stage].lift_index = next_lift_index - 1;
        grid_[stage].type = GridType::Intermediate;
        if (numerics::isApprox(ti+dt, next_impulse_time, eps)) {
          ti += dt;
          grid_[stage].dt = ti + dt - next_impulse_time;
        }
      }
    }
    if (has_next_lift) {
      const double next_lift_time = contact_sequence->liftTime(next_lift_index);
      if (next_lift_time <= ti+dt+eps && (next_lift_time+margin < t+T_)) {
        grid_[stage].dt = next_lift_time - ti;
        ++stage;
        ++next_lift_index;
        grid_[stage].t = next_lift_time;
        grid_[stage].dt = std::min(ti+dt, t+T_) - next_lift_time;
        grid_[stage].time_stage = stage;
        grid_[stage].contact_phase = next_impulse_index + next_lift_index;
        grid_[stage].impulse_index = next_impulse_index - 1;
        grid_[stage].lift_index = next_lift_index - 1;
        grid_[stage].type = GridType::Lift;
        if (numerics::isApprox(ti+dt, next_lift_time, eps)) {
          ti += dt;
          grid_[stage].dt = ti + dt - next_lift_time;
        }
      }
    }
    ++stage;
    ti += dt;
  }
  grid_[stage].t  = t+T_;
  grid_[stage].dt = 0;
  grid_[stage].time_stage = stage;
  grid_[stage].contact_phase = next_impulse_index + next_lift_index;
  grid_[stage].impulse_index = next_impulse_index - 1;
  grid_[stage].lift_index = next_lift_index - 1;
  grid_[stage].type = GridType::Terminal;
  N_grids_ = stage;

  // set dt_next
  for (int i=0; i<N_grids_; ++i) {
    grid_[i].dt_next = grid_[i+1].dt;
  }
  grid_[N_grids_].dt_next = 0.0;

  // set switching_constraint flag
  for (int i=0; i<N_grids_-1; ++i) {
    grid_[i].switching_constraint = (grid_[i+2].type == GridType::Impulse);
  }
  grid_[N_grids_-1].switching_constraint = false;
  grid_[N_grids_].switching_constraint = false;

  // set t0 and disable sto
  for (int i=0; i<=N_grids_; ++i) {
    grid_[i].t0 = t;
    grid_[i].sto = false;
    grid_[i].sto_next = false;
    grid_[i].grid_count_in_phase = 1;
    grid_[i].N_phase = 1;
  }

  // count grids
  int grid_count_in_phase = 0;
  int phase_start_stage = 0;
  for (int i=0; i<N_grids_; ++i) {
    if (grid_[i].type == GridType::Impulse) {
      for (int j=phase_start_stage; j<i; ++j) {
        grid_[j].N_phase = grid_count_in_phase;
      }
      grid_[i].grid_count_in_phase = 0;
      grid_[i].N_phase = 0;
      ++i;
      grid_count_in_phase = 0;
      phase_start_stage = i;
    }
    else if (grid_[i].type == GridType::Lift) {
      for (int j=phase_start_stage; j<i; ++j) {
        grid_[j].N_phase = grid_count_in_phase;
      }
      grid_count_in_phase = 0;
      phase_start_stage = i;
    }
    grid_[i].grid_count_in_phase = grid_count_in_phase;
    ++grid_count_in_phase;
  }
  for (int j=phase_start_stage; j<N_grids_; ++j) {
    grid_[j].N_phase = grid_count_in_phase;
  }
  grid_[N_grids_].grid_count_in_phase = 0;
  grid_[N_grids_].N_phase = 0;
}


void TimeDiscretization::discretizePhase(
    const std::shared_ptr<ContactSequence>& contact_sequence, const double t) {
  int prev_event_stage = 0;
  double prev_event_time = t;
  for (int i=0; i<N_grids_; ++i) {
    if (grid_[i].type == GridType::Impulse) {
      const double event_time = contact_sequence->impulseTime(grid_[i+1].impulse_index);
      const double dt = (event_time - prev_event_time) / grid_[i-1].N_phase;
      for (int j=prev_event_stage; j<=i-1; ++j) {
        grid_[j].t = prev_event_time + (j-prev_event_stage) * dt;
        grid_[j].dt = dt;
      }
      grid_[i].t = event_time; 
      grid_[i].dt = 0.0; 
      prev_event_time = event_time;
      prev_event_stage = i+1;
      ++i;
    }
    else if (grid_[i+1].type == GridType::Lift) {
      const double event_time = contact_sequence->liftTime(grid_[i+1].lift_index);
      const double dt = (event_time - prev_event_time) / grid_[i].N_phase;
      for (int j=prev_event_stage; j<=i; ++j) {
        grid_[j].t = prev_event_time + (j-prev_event_stage) * dt;
        grid_[j].dt = dt;
      }
      prev_event_time = event_time;
      prev_event_stage = i+1;
    }
    else if (grid_[i+1].type == GridType::Terminal) {
      const double dt = (t + T_ - prev_event_time) / grid_[i].N_phase;
      for (int j=prev_event_stage; j<=i; ++j) {
        grid_[j].t = prev_event_time + (j-prev_event_stage) * dt;
        grid_[j].dt = dt;
      }
    }
  }
  grid_[N_grids_].t = t + T_;
  grid_[N_grids_].dt = 0;

  // set dt_next
  for (int i=0; i<N_grids_; ++i) {
    grid_[i].dt_next = grid_[i+1].dt;
  }
  grid_[N_grids_].dt_next = 0.0;

  // set sto
  std::vector<bool> sto_event;
  sto_event.reserve(grid_[N_grids_].contact_phase-grid_[N_grids_].contact_phase);
  for (int i=0; i<N_grids_; ++i) {
    if (grid_[i].type == GridType::Impulse) {
      sto_event.push_back(
          contact_sequence->isSTOEnabledImpulse(grid_[i+1].impulse_index));
    }
    else if (grid_[i].type == GridType::Lift) {
      sto_event.push_back(
          contact_sequence->isSTOEnabledLift(grid_[i+1].lift_index));
    }
  }
  if (sto_event.empty()) return;

  std::vector<bool> sto_phase;
  sto_phase.reserve(grid_[N_grids_].contact_phase-grid_[N_grids_].contact_phase+2);
  sto_phase.push_back(sto_event.front());
  for (int i=1; i<sto_event.size(); ++i) {
    sto_phase.push_back(sto_event[i-1] || sto_event[i]);
  }
  sto_phase.push_back(sto_event.back());
  sto_phase.push_back(false);

  for (int i=0; i<N_grids_; ++i) {
    const int phase = grid_[i].contact_phase - grid_[0].contact_phase;
    grid_[i].sto = sto_phase[phase];
    grid_[i].sto_next = sto_phase[phase+1];
  }
  grid_[N_grids_].sto = false;
  grid_[N_grids_].sto_next = false;
}


void TimeDiscretization::disp(std::ostream& os) const {
  os << "Time discretization of optimal control problem (OCP):" << "\n";
  os << "  T: " << T_ << "\n";
  os << "  N_ideal: " << N_ideal() << "\n";
  os << "  N: " << N() << "\n";
  os << "  N_impulse: " << N_impulse() << "\n";
  os << "  N_lift: " << N_lift() << "\n";
  os << "  N_all: " << (N()+1+2*N_impulse()+N_lift()) << "\n";
  os << "  No. of discrete events: " << numDiscreteEvents() << "\n";
  os << "  No. of contact phases: " << numContactPhases() << "\n";
  for (int i=0; i<numContactPhases(); ++i) {
    os << "    No. of grids at contact phase " << i << ": " << N_phase(i) << "\n";
  }
  os << "  isFormulationTractable: " << std::boolalpha 
     << isFormulationTractable() << "\n";;
  os << "  isSwitchingTimeConsistent: " << std::boolalpha
     << isSwitchingTimeConsistent() << "\n";
  os << " -----------------------------------------------------------------------" << "\n";
  os << "  grid point | grid count |      t |     dt | phase |  sto  | sto_next |" << "\n";
  os << " -----------------------------------------------------------------------" << "\n";
  for (int i=0; i<N(); ++i) {
    os << "  stage: " << std::right << std::setw(3) << i << " | "
       << "       " << std::right << std::setw(3) << gridInfo(i).grid_count_in_phase << " | " 
       << std::fixed << std::setprecision(4) << gridInfo(i).t << " | " << gridInfo(i).dt
       << " |   " << std::setw(3) << contactPhase(i)
       << " | " << std::setw(5) << std::boolalpha << isSTOEnabledPhase(contactPhase(i)) 
       << " |   " << std::setw(5) << std::boolalpha << isSTOEnabledNextPhase(contactPhase(i)) 
       << "  |" << "\n";
    if (isTimeStageBeforeImpulse(i)) {
      const int impulse_index = impulseIndexAfterTimeStage(i);
      os << "    aux: " << std::right << std::setw(3) << impulse_index << " | "
         << std::right << std::setw(3) << "         0" << " | " 
         << std::fixed << std::setprecision(4) << gridInfoImpulse(impulse_index).t
         << " | " << gridInfoImpulse(impulse_index).dt
         << " |   " << std::setw(3) << contactPhaseAfterImpulse(impulse_index)
         << " | " << std::setw(5) << std::boolalpha 
         << isSTOEnabledPhase(contactPhaseAfterImpulse(impulse_index)) 
         << " |   " << std::setw(5) << std::boolalpha 
         << isSTOEnabledNextPhase(contactPhaseAfterImpulse(impulse_index)) 
         << "  |" << "\n";
    }
    else if (isTimeStageBeforeLift(i)) {
      const int lift_index = liftIndexAfterTimeStage(i);
      os << "   lift: " << std::right << std::setw(3) << lift_index << " | "
         << std::right << std::setw(3) << "         0" << " | " 
         << std::fixed << std::setprecision(4) << gridInfoLift(lift_index).t
         << " | " << gridInfoLift(lift_index).dt
         << " |   " << std::setw(3) << contactPhaseAfterLift(lift_index)
         << " | " << std::setw(5) << std::boolalpha
         << isSTOEnabledPhase(contactPhaseAfterLift(lift_index)) 
         << " |   " << std::setw(5) << std::boolalpha
         << isSTOEnabledNextPhase(contactPhaseAfterLift(lift_index))
         << "  |" << "\n";
    }
  }
  os << "  stage: " << std::right << std::setw(3) << N() << " | " 
     << "       " << std::right << std::setw(3) << gridInfo(N()).grid_count_in_phase << " | " 
     << std::fixed << std::setprecision(4) << gridInfo(N()).t
     << " |        |       |       |          |" << "\n"; 
  os << " -----------------------------------------------------------------------" << std::flush;
}


std::ostream& operator<<(std::ostream& os, 
                         const TimeDiscretization& time_discretization) {
  time_discretization.disp(os);
  return os;
}

} // namespace robotoc 