#include "urdf_factory.hpp"


namespace robotoc {
namespace testhelper {

std::string RobotManipulatorURDF() {
  return "../urdf/iiwa14/iiwa14.urdf";
}


std::string QuadrupedURDF() {
  return "../urdf/anymal/anymal.urdf";
}


std::string HumanoidURDF() {
  return "../urdf/icub/icub.urdf";
}

} // namespace testhelper
} // namespace robotoc