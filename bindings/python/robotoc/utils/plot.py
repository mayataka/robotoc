import numpy as np
import sys
import os


class PlotConvergence:
    def __init__(self):
        self.wspace = 0.45
        self.hspace = 0.0
        self.figsize = 3, 2
        self.figsize_sto = 6, 2
        self.legend_bbox_to_anchor = (0.075, 1.0)
        self.legend_ncols = 2
        self.ylim = None
        self.xlim = None
        self.linestyle = ['solid', 'dashed', 'dashdot', 
                          (0, (3, 1, 1, 1, 1, 1)),
                          (0, (3, 1, 1, 1, 1))]


    def plot(self, kkt_data, ts_data=None, fig_name=None, save_dir='./'):
        import matplotlib.pyplot as plt
        import seaborn 
        seaborn.set()
        seaborn.set_style("whitegrid")
        seaborn.set_style("ticks")
        seaborn.set_palette("deep")
        seaborn.set_context("paper")
        plt.rc('mathtext', 
            **{'rm':'serif', 
            'it':'serif:itelic', 
            'bf':'serif:bold', 
            'fontset':'cm'}
        )
        plt.rcParams['xtick.direction'] = 'in'
        plt.rcParams['ytick.direction'] = 'in'
        plt.rcParams['pdf.fonttype'] = 42
        plt.rcParams['ps.fonttype'] = 42
        # Set the width of lines and axes.
        plt.rcParams['lines.linewidth'] = 1.5
        plt.rcParams['axes.linewidth'] = 0.5
        if ts_data is not None:
            plt.rcParams['figure.figsize'] = self.figsize_sto
        else:
            plt.rcParams['figure.figsize'] = self.figsize

        fig = plt.figure()
        if ts_data is not None:
            ax_ts = fig.add_subplot(1, 2, 1)
            ax_kkt = fig.add_subplot(1, 2, 2)
            for i, ts in enumerate(np.array(ts_data).T):
                ts_str = r"$t_{}$".format(i+1)
                linestyle = self.linestyle[i%len(self.linestyle)]
                ax_ts.plot(ts, label=ts_str+r"${\rm [s]}$", linestyle=linestyle)
            if self.ylim is not None:
                ax_ts.set_ylim(self.ylim)
            ax_ts.set_xlabel("No. of Iterations")
            ax_ts.set_ylabel(r"$t_k {\rm [s]}$")
            ax_ts.legend(
                edgecolor='black', loc='upper left', 
                bbox_to_anchor=self.legend_bbox_to_anchor, 
                ncol=self.legend_ncols
            )
        else:
            ax_kkt = fig.add_subplot(1, 1, 1)

        ax_kkt.plot(np.log10(kkt_data))
        ax_kkt.set_xlabel("No. of Iterations")
        ax_kkt.set_ylabel(r"$\log_{10} \|$ (KKT residual) $\|_2$")

        os.makedirs(os.path.abspath(save_dir), exist_ok=True)
        if fig_name is not None:
            plt.subplots_adjust(wspace=self.wspace, hspace=self.hspace)
            plt.savefig(
                os.path.join(save_dir, fig_name+'.pdf'),
                bbox_inches="tight", 
                pad_inches=0.1
            )
        else:
            plt.show()


class PlotContactForce:
    def __init__(self, mu=None):
        self.mu = mu
        self.eps = sys.float_info.epsilon
        self.wspace = 0.3
        self.hspace = 0.45
        self.figsize = 6, 4
        self.legend_bbox_to_anchor = (-0.2, 3.0)
        self.ylim = [-200, 200]

    def _decompose_f_data(self, f_data):
        length = len(f_data) 
        fx = []
        fy = []
        fz = []
        mfz = []
        for i in range(length):
            fx.append(f_data[i][0])
            fy.append(f_data[i][1])
            fz.append(f_data[i][2])
            mfz.append(-1*f_data[i][2])
        return fx, fy, fz, mfz

    def _detect_phases(self, fx, fy, fz, mfz):
        assert len(fx) == len(fy)
        assert len(fx) == len(fz)
        assert len(fx) == len(mfz)
        N = len(fx) 
        f_norm = [fz[i]*fz[i] for i in range(N)]
        N_phase = []
        active_phase = []
        active = f_norm[0] > self.eps
        ngrids = 0
        for i in range(N):
            if active:
                if f_norm[i] <= self.eps:
                    N_phase.append(ngrids)
                    active_phase.append(active)
                    active = False
                    ngrids = 0
                else:
                    ngrids += 1
            else:
                if f_norm[i] > self.eps:
                    N_phase.append(ngrids)
                    active_phase.append(active)
                    active = True
                    ngrids = 0
                else:
                    ngrids += 1
        N_phase.append(ngrids)
        active_phase.append(active)
        return N_phase, active_phase

    def _plot_f_leg(self, ax, fx, fy, fz, mfz, t, title, enable_xlabel):
        import matplotlib.pyplot as plt
        cmap = plt.get_cmap("tab10")
        N_phase, active_phase = self._detect_phases(fx, fy, fz, mfz)
        for phase in range(len(N_phase)):
            if active_phase[phase]:
                N_begin = 0
                if phase > 0:
                    N_begin = sum([N_phase[i] for i in range(phase)]) + 1
                N_end = N_begin + N_phase[phase] 
                if phase == len(N_phase)-1:
                    ax.plot(
                        t[N_begin:N_end], 
                        fx[N_begin:N_end], 
                        color=cmap(0), 
                        linewidth=1.5, 
                        linestyle='-'
                    )
                    ax.plot(
                        t[N_begin:N_end], 
                        fy[N_begin:N_end], 
                        color=cmap(2), 
                        linewidth=1.5, 
                        linestyle='--'
                    )
                else:
                    ax.plot(
                        t[N_begin:N_end], 
                        fx[N_begin:N_end], 
                        color=cmap(0), 
                        linewidth=1.5, 
                        linestyle='-', 
                        label=r'$f_x {\rm [N]}$'
                    )
                    ax.plot(
                        t[N_begin:N_end], 
                        fy[N_begin:N_end], 
                        color=cmap(2), 
                        linewidth=1.5, 
                        linestyle='--', 
                        label=r'$f_y {\rm [N]}$'
                    )
                ax.set_ylim(self.ylim)
                if self.mu is not None:
                    ylims = ax.get_ylim()
                    llim = np.full(len(fz), ylims[0]) 
                    ulim = np.full(len(fz), ylims[1]) 
                    mufz = [self.mu / np.sqrt(2) * f for f in fz]
                    mumfz = [self.mu / np.sqrt(2) * mf for mf in mfz]
                    ax.fill_between(
                        t[N_begin:N_end], 
                        mufz[N_begin:N_end], 
                        ulim[N_begin:N_end], 
                        color='gray', 
                        edgecolor='black', 
                        alpha=0.25, 
                        hatch="///////"
                    )
                    ax.fill_between(
                        t[N_begin:N_end], 
                        mumfz[N_begin:N_end], 
                        llim[N_begin:N_end], 
                        color='gray', 
                        edgecolor='black', 
                        alpha=0.25, 
                        hatch="///////"
                    )
                if phase == len(N_phase)-1:
                    ax.plot(
                        t[N_begin:N_end], 
                        fz[N_begin:N_end], 
                        color=cmap(1), 
                        linewidth=1.5, 
                        linestyle='-.'
                    )
                else:
                    ax.plot(
                        t[N_begin:N_end], 
                        fz[N_begin:N_end], 
                        color=cmap(1), 
                        linewidth=1.5, 
                        linestyle='-.', 
                        label=r'$f_z {\rm [N]}$'
                    )
        if enable_xlabel:
            ax.set_xlabel(r'${\rm Time \; [s]}$')
        ax.set_title(title)
        ax.set_xlim([t[0], t[-1]])


    def plot(self, f_data, t, fig_name=None, save_dir='./'):
        import matplotlib.pyplot as plt
        import seaborn 
        assert len(t) >= len(f_data)
        seaborn.set()
        seaborn.set_style('whitegrid')
        seaborn.set_style("ticks")
        seaborn.set_palette("deep") 
        seaborn.set_context("paper")
        plt.rc('mathtext', 
            **{'rm':'serif', 
            'it':'serif:itelic', 
            'bf':'serif:bold', 
            'fontset':'cm'}
        )
        plt.rcParams['xtick.direction'] = 'in'
        plt.rcParams['ytick.direction'] = 'in'
        plt.rcParams['pdf.fonttype'] = 42
        plt.rcParams['ps.fonttype'] = 42
        # Set the width of lines and axes.
        plt.rcParams['lines.linewidth'] = 1
        plt.rcParams['axes.linewidth'] = 0.5
        plt.rcParams['figure.figsize'] = self.figsize

        plt.rcParams['xtick.labelsize'] = 8
        plt.rcParams['ytick.labelsize'] = 8

        fig = plt.figure()
        ax_LF = fig.add_subplot(2, 2, 1)
        ax_LH = fig.add_subplot(2, 2, 2)
        ax_RF = fig.add_subplot(2, 2, 3)
        ax_RH = fig.add_subplot(2, 2, 4)

        plt.subplots_adjust(wspace=self.wspace, hspace=self.hspace) 

        f_LF = []
        f_LH = []
        f_RF = []
        f_RH = []
        N = len(f_data) 
        for i in range(N):
            f_LF.append(f_data[i][0:3])
            f_LH.append(f_data[i][3:6])
            f_RF.append(f_data[i][6:9])
            f_RH.append(f_data[i][9:12])

        fx_LF, fy_LF, fz_LF, mfz_LF = self._decompose_f_data(f_LF)
        fx_LH, fy_LH, fz_LH, mfz_LH = self._decompose_f_data(f_LH)
        fx_RF, fy_RF, fz_RF, mfz_RF = self._decompose_f_data(f_RF)
        fx_RH, fy_RH, fz_RH, mfz_RH = self._decompose_f_data(f_RH)

        self._plot_f_leg(ax_LF, fx_LF, fy_LF, fz_LF, mfz_LF, t, 'Left-front leg', False)
        self._plot_f_leg(ax_LH, fx_LH, fy_LH, fz_LH, mfz_LH, t, 'Left-hip leg', False)
        self._plot_f_leg(ax_RF, fx_RF, fy_RF, fz_RF, mfz_RF, t, 'Right-front leg', True)
        self._plot_f_leg(ax_RH, fx_RH, fy_RH, fz_RH, mfz_RH, t, 'Right-hip leg', True)

        plt.legend(
            edgecolor='black', 
            loc='upper center', 
            bbox_to_anchor=self.legend_bbox_to_anchor, 
            ncol=4
        )

        os.makedirs(os.path.abspath(save_dir), exist_ok=True)
        if fig_name is not None:
            plt.savefig(
                os.path.join(save_dir, fig_name+'.pdf'),
                bbox_inches="tight", 
                pad_inches=0.1
            )
        else:
            plt.show()


class PlotLeggedMPC:
    def __init__(self):
        self.wspace = 0.45
        self.hspace = 0.1
        self.figsize = 6, 2
        self.legend_bbox_to_anchor = (0.515, 1.0)
        self.legend_ncols = 2
        self.ylim = None
        self.xlim = None

    def plot(self, t_data, vcom_data, wcom_data, 
             vcom_command=None, yaw_rate_command=None, 
             fig_name=None, save_dir='./'):
        import matplotlib.pyplot as plt
        import seaborn 
        seaborn.set()
        seaborn.set_style("whitegrid")
        seaborn.set_style("ticks")
        seaborn.set_palette("deep")
        seaborn.set_context("paper")
        plt.rc('mathtext', 
            **{'rm':'serif', 
            'it':'serif:itelic', 
            'bf':'serif:bold', 
            'fontset':'cm'}
        )
        plt.rcParams['xtick.direction'] = 'in'
        plt.rcParams['ytick.direction'] = 'in'
        plt.rcParams['pdf.fonttype'] = 42
        plt.rcParams['ps.fonttype'] = 42
        # Set the width of lines and axes.
        plt.rcParams['lines.linewidth'] = 1.5
        plt.rcParams['axes.linewidth'] = 0.5

        cmap = plt.get_cmap("tab10")
        fig = plt.figure()
        ax_v = fig.add_subplot(2, 1, 1)
        ax_w = fig.add_subplot(2, 1, 2)
        ax_v.plot(t_data, np.array(vcom_data).T[0], label=r"$v_x {\rm [m/s]}$", 
                  linestyle='-', color=cmap(0))
        ax_v.plot(t_data, np.array(vcom_data).T[1], label=r"$v_y {\rm [m/s]}$", 
                  linestyle='-', color=cmap(1))
        ax_v.plot(t_data, np.array(vcom_data).T[2], label=r"$v_z {\rm [m/s]}$", 
                  linestyle='-', color=cmap(2))
        if vcom_command is not None:
            ax_v.plot(t_data, np.array(vcom_command).T[0], 
                      label=r"$v_{x, {\rm cmd}} {\rm [m/s]}$", 
                      linestyle='--', color=cmap(3))
            ax_v.plot(t_data, np.array(vcom_command).T[1], 
                      label=r"$v_{y, {\rm cmd}} {\rm [m/s]}$", 
                      linestyle='--', color=cmap(4))
        ax_w.plot(t_data, wcom_data)
        ax_w.plot(t_data, np.array(wcom_data).T[0], label=r"$w_x {\rm [rad/s]}$", 
                  linestyle='-', color=cmap(0))
        ax_w.plot(t_data, np.array(wcom_data).T[1], label=r"$w_y {\rm [rad/s]}$", 
                  linestyle='-', color=cmap(1))
        ax_w.plot(t_data, np.array(wcom_data).T[2], label=r"$w_z {\rm [rad/s]}$", 
                  linestyle='-', color=cmap(2))
        if yaw_rate_command is not None:
            ax_w.plot(t_data, yaw_rate_command, 
                      label=r"$w_{z, {\rm cmd}} {\rm [rad/s]}$", 
                      linestyle='--', color=cmap(5))
        ax_v.set_xticklabels([])
        ax_v.set_ylabel("Linear CoM velocity")
        ax_w.set_ylabel("Angular CoM velocity")
        ax_w.set_xlabel(r"${\rm Time \; [s]}$")
        ax_v.legend(
            edgecolor='black', loc='upper left', 
            bbox_to_anchor=self.legend_bbox_to_anchor, 
            ncol=self.legend_ncols
        )
        ax_w.legend(
            edgecolor='black', loc='upper left', 
            bbox_to_anchor=self.legend_bbox_to_anchor, 
            ncol=self.legend_ncols
        )

        os.makedirs(os.path.abspath(save_dir), exist_ok=True)
        if fig_name is not None:
            plt.subplots_adjust(wspace=self.wspace, hspace=self.hspace)
            plt.savefig(
                os.path.join(save_dir, fig_name+'.pdf'),
                bbox_inches="tight", 
                pad_inches=0.1
            )
        else:
            plt.show()