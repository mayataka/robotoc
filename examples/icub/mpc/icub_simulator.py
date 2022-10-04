import robotoc_sim


class iCubSimulator(robotoc_sim.LeggedSimulator):
    def __init__(self, urdf_path, time_step):
        super().__init__(urdf_path, time_step)

    @classmethod
    def get_joint_id_map(self):
        return [1, 2, 4, 5, 6, 7, 12, 13, 15, 16, 17, 18]