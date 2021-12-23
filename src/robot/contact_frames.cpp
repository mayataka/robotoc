#include "robotoc/robot/contact_frames.hpp"


namespace robotoc {

ContactFrames::ContactFrames(const std::vector<int>& _point_contact_frames, 
                             const std::vector<int>& _surface_contact_frames)
  : point_contact_frames(_point_contact_frames),
    surface_contact_frames(_surface_contact_frames) {
}


ContactFrames::ContactFrames()
  : point_contact_frames(),
    surface_contact_frames() {
}


ContactFrames::~ContactFrames() {
}

} // namespace robotoc 