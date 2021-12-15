import time


def cpu_time(ocp_solver, t, q, v, num_iteration):
    start_clock = time.time()
    for i in range(num_iteration):
        ocp_solver.update_solution(t, q, v)
    end_clock = time.time()
    print('---------- OCP benchmark : CPU time ----------')
    print('total CPU time: {:.6g}'.format(1e03*(end_clock-start_clock)) + '[ms]')
    print('CPU time per update: {:.6g}'.format(1e03*(end_clock-start_clock)/num_iteration) + '[ms]')
    print('-----------------------------------')