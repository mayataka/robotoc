import time


def cpu_time(ocp_solver, t, q, v, num_iteration, line_search=False):
    start_clock = time.time()
    for i in range(num_iteration):
        ocp_solver.update_solution(t, q, v, line_search)
    end_clock = time.time()
    print('---------- OCP benchmark : CPU time ----------')
    print('total CPU time: {:.6g}'.format(1e03*(end_clock-start_clock)) + '[ms]')
    print('CPU time per update: {:.6g}'.format(1e03*(end_clock-start_clock)/num_iteration) + '[ms]')
    print('-----------------------------------')


def convergence(ocp_solver, t, q, v, num_iteration, line_search=False):
    print('---------- OCP benchmark : Convergence ----------')
    ocp_solver.compute_KKT_residual(t, q, v)
    print('Initial KKT error = {:.6g}'.format(ocp_solver.KKT_error()))
    for i in range(num_iteration):
        ocp_solver.update_solution(t, q, v, line_search)
        ocp_solver.compute_KKT_residual(t, q, v)
        print('KKT error after iteration ' + str(i+1) + ' = {:.6g}'.format(ocp_solver.KKT_error()))
    print('-----------------------------------')