import pybullet
import numpy as np
import os


def get_control_input(control_policy, q: np.ndarray, v: np.ndarray):
    nJ = control_policy.tauJ.shape[0]
    qJ = q[-nJ:]
    dqJ = v[-nJ:]
    return control_policy.tauJ - control_policy.Kp @ (control_policy.qJ-qJ) \
                               - control_policy.Kd @ (control_policy.dqJ-dqJ)

class MPCSimulation(object):
    def __init__(self, simulator):
        self.simulator = simulator

    def run(self, mpc, t0: float, q0: np.ndarray, simulation_time: float, 
            feedback_policy: bool=False, feedback_delay: bool=False, 
            simulation_steps_per_mpc_update: int=1, verbose: bool=False, 
            log: bool=False, record: bool=False, name: str='mpc_sim'):
        assert simulation_steps_per_mpc_update >= 1
        self.simulator.init_simulation(t0, q0)
        self.name = name
        if log:
            log_dir = os.path.join(os.getcwd(), name+"_log")
            self.log_dir = log_dir
            os.makedirs(log_dir, exist_ok=True)
            q_log = open(os.path.join(log_dir, "q.log"), mode='w')
            v_log = open(os.path.join(log_dir, "v.log"), mode='w')
            u_log = open(os.path.join(log_dir, "u.log"), mode='w')
            t_log = open(os.path.join(log_dir, "t.log"), mode='w')
            kkt_log = open(os.path.join(log_dir, "kkt.log"), mode='w')
        if record:
            pybullet.startStateLogging(pybullet.STATE_LOGGING_VIDEO_MP4, name+".mp4")

        use_feedback_policy = (simulation_steps_per_mpc_update >= 2) and feedback_policy
        inner_loop_count = simulation_steps_per_mpc_update - 1
        dt = self.simulator.time_step
        while self.simulator.get_time() < t0 + simulation_time:
            t = self.simulator.get_time()
            q, v = self.simulator.get_state()
            if verbose:
                print('t = {:.6g}:'.format(t))
            if feedback_delay:
                if use_feedback_policy:
                    u = get_control_input(mpc.get_control_policy(t), q, v) 
                else:
                    u = mpc.get_initial_control_input().copy()
            if inner_loop_count == 0:
                mpc.update_solution(t, dt, q, v)
                inner_loop_count = simulation_steps_per_mpc_update - 1
            else:
                inner_loop_count -= 1
            kkt_error = mpc.KKT_error(t, q, v) 
            if verbose:
                print('KKT error = {:.6g}'.format(kkt_error))
                print('')
            if not feedback_delay:
                if use_feedback_policy:
                    u = get_control_input(mpc.get_control_policy(t), q, v) 
                else:
                    u = mpc.get_initial_control_input().copy()
            self.simulator.step_simulation(u)
            if log:
                np.savetxt(q_log, [q])
                np.savetxt(v_log, [v])
                np.savetxt(u_log, [u])
                np.savetxt(t_log, np.array([t]))
                np.savetxt(kkt_log, np.array([kkt_error]))

        if log:
            q_log.close()
            v_log.close()
            u_log.close()
            t_log.close()
            kkt_log.close()
            self.q_log = os.path.abspath(os.path.join(log_dir, "q.log"))
            self.v_log = os.path.abspath(os.path.join(log_dir, "v.log"))
            self.u_log = os.path.abspath(os.path.join(log_dir, "u.log"))
            self.t_log = os.path.abspath(os.path.join(log_dir, "t.log"))
            self.kkt_log = os.path.abspath(os.path.join(log_dir, "kkt.log"))
            self.log_dir = os.path.abspath(log_dir)
            print('Logs are saved at ' + self.log_dir)

        if record:
            self.simulator.disconnect()
