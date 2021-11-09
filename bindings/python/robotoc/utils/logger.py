import os
import numpy as np
import robotoc


class Logger:
    def __init__(self, vars, log_name, root_dir='./'):
        self.vars = vars
        self.log_dir = os.path.abspath(root_dir+log_name+'_log') 
        self.logs = [os.path.join(self.log_dir, var+'.log') for var in vars]
        os.makedirs(self.log_dir, exist_ok=True)
        for log in self.logs:
            with open(log, mode='w') as f:
                f.write('')
                f.close()

    def take_log(self, solver, precision='%.18e', delimiter=', '):
        for var, log in zip(self.vars, self.logs):
            if var == 'KKT':
                with open(log, 'a') as logf:
                    np.savetxt(logf, np.array([solver.KKT_error()]), delimiter=delimiter)
            elif var == 'ts' and isinstance(solver, robotoc.OCPSolver):
                ts = solver.get_solution('ts')
                with open(log, 'a') as logf:
                    np.savetxt(logf, ts, fmt=precision, delimiter=delimiter)
            elif var == 'f' and isinstance(solver, robotoc.OCPSolver):
                fs = solver.get_solution('f', 'WORLD')
                with open(log, 'a') as logf:
                    np.savetxt(logf, fs, fmt=precision, delimiter=delimiter)
            else:
                sols = solver.get_solution(var)
                with open(log, 'a') as logf:
                    np.savetxt(logf, sols, fmt=precision, delimiter=delimiter)

    def get_data(self, var):
        log_file = os.path.join(self.log_dir, var+'.log')  
        return np.genfromtxt(log_file, delimiter=',')