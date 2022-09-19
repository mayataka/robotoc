import robotoc_sim


class A1Simulator(robotoc_sim.LeggedSimulator):
    def __init__(self, urdf_path, time_step):
        super().__init__(urdf_path, time_step)

    @classmethod
    def get_joint_id_map(self):
        return [7, 9, 10, 2, 4, 5, 17, 19, 20, 12, 14, 15]