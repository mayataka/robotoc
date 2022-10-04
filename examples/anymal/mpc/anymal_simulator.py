import robotoc_sim


class ANYmalSimulator(robotoc_sim.LeggedSimulator):
    def __init__(self, urdf_path, time_step):
        super().__init__(urdf_path, time_step)

    @classmethod
    def get_joint_id_map(self):
        return [1, 2, 3, 11, 12, 13, 6, 7, 8, 16, 17, 18]