import robotoc_sim


class ANYmalCSimulator(robotoc_sim.LeggedSimulator):
    def __init__(self, urdf_path, time_step):
        super().__init__(urdf_path, time_step)

    @classmethod
    def get_joint_id_map(self):
        return [37, 40, 43, 57, 60, 63, 47, 50, 53, 67, 70, 73]