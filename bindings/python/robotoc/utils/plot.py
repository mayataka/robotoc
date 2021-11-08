import matplotlib.pyplot as plt
import seaborn 
import numpy as np



class PlotContactForce:
    def __init__(self, mu=None):
        self.mu = mu

    def _decompose_f(self, f_data):
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

    def _plot_f_leg(self, ax, f, N, t, title, enable_xlabel, ylim=None, fill_feas=False):
        assert len(N) == len(t)
        cmap = plt.get_cmap("tab10")
        for Ni, ti in zip(N, t):
            ax.plot(t1, fx[:Ni], color=cmap(0), linewidth=1.5, linestyle='-')
            ax.plot(t2, fx[-N3:], color=cmap(0), linewidth=1.5, label=r'$f_x {\rm [N]}$', linestyle='-')
            ax.plot(t1, fy[:Ni], color=cmap(2), linewidth=1.5, linestyle='--')
            ax.plot(t2, fy[-N3:], color=cmap(2), linewidth=1.5, label=r'$f_y {\rm [N]}$', linestyle='--')
        ax.plot(np.linspace(0, dt1*N1+dt2*N2+dt3*N3, N), np.zeros(N), color='black', alpha=0.25, linestyle=':', linewidth=1.2)
        if enable_xlabel:
            ax.set_xlabel(r'${\rm Time \; [s]}$')
        ax.set_title(title)
        if ylim is not None:
            ax.set_ylim(ylim)
        else:
            ax.set_ylim([-400, 400])
        if fill_feas:
            ylims = ax.get_ylim()
            llim1 = np.full(N1, ylims[0]) 
            ulim1 = np.full(N1, ylims[1]) 
            llim3 = np.full(N3, ylims[0]) 
            ulim3 = np.full(N3, ylims[1]) 
            ax.fill_between(t1, fz[:N1], ulim1, color='gray', edgecolor='black', alpha=0.25, hatch="///////")
            ax.fill_between(t1, mfz[:N1], llim1, color='gray', edgecolor='black', alpha=0.25, hatch="///////")
            ax.fill_between(t2, fz[-N3:], ulim3, color='gray', edgecolor='black', alpha=0.25, hatch="///////")
            ax.fill_between(t2, mfz[-N3:], llim3, color='gray', edgecolor='black', alpha=0.25, hatch="///////")
            # ax.fill_between(t1, fz[:20], ulim, color='gray', alpha=0.3, hatch="//////")
            # ax.fill_between(t1, mfz[:20], llim, color='gray', alpha=0.3, hatch="//////")
            # ax.fill_between(t2, fz[20:], ulim, color='gray', alpha=0.3, hatch="//////")
            # ax.fill_between(t2, mfz[20:], llim, color='gray', alpha=0.3, hatch="//////")
            mufz = [np.sqrt(2) * f / self.mu for f in fz]
            ax.plot(t1, mufz[:N1], color=cmap(1), linewidth=1.5, linestyle='-.')
            ax.plot(t2, mufz[-N3:], color=cmap(1), linewidth=1.5, label=r'$f_z {\rm [N]}$', linestyle='-.')
        else:
            ax.plot(t1, fz[:N1], color=cmap(1), alpha=0.5, linestyle='-.')
            ax.plot(t2, fz[-N3:], color=cmap(1), alpha=0.5, label=r'$\frac{\mu}{\sqrt{2}} f_z {\rm [N]}$', linestyle='-.')
            ax.plot(t1, mfz[:N1], color=cmap(3), alpha=0.5, linestyle='-.')
            ax.plot(t2, mfz[-N3:], color=cmap(3), alpha=0.5, label=r'$- \frac{\mu}{\sqrt{2}} f_z {\rm [N]}$', linestyle='-.')

    def plot_f(self, f_data, N1, N2, N3, dt1, dt2, dt3, fig_name):
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
        plt.rcParams['figure.figsize'] = 6, 4

        plt.rcParams['xtick.labelsize'] = 8
        plt.rcParams['ytick.labelsize'] = 8

        fig = plt.figure()
        ax_LF = fig.add_subplot(2, 2, 1)
        ax_LH = fig.add_subplot(2, 2, 2)
        ax_RF = fig.add_subplot(2, 2, 3)
        ax_RH = fig.add_subplot(2, 2, 4)

        plt.subplots_adjust(wspace=0.3, hspace=0.45) 

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

        fx_LF, fy_LF, fz_LF, mfz_LF = self._decompose_f(f_LF)
        fx_LH, fy_LH, fz_LH, mfz_LH = self._decompose_f(f_LH)
        fx_RF, fy_RF, fz_RF, mfz_RF = self._decompose_f(f_RF)
        fx_RH, fy_RH, fz_RH, mfz_RH = self._decompose_f(f_RH)

        self.plot_f_leg(ax_LF, fx_LF, fy_LF, fz_LF, mfz_LF, N1, N2, N3, dt1, dt2, dt3, 'Left-front leg', False, ylim=[-250, 250], fill_feas=True)
        self.plot_f_leg(ax_LH, fx_LH, fy_LH, fz_LH, mfz_LH, N1, N2, N3, dt1, dt2, dt3, 'Left-hip leg', False, ylim=[-250, 250], fill_feas=True)
        self.plot_f_leg(ax_RF, fx_RF, fy_RF, fz_RF, mfz_RF, N1, N2, N3, dt1, dt2, dt3, 'Right-front leg', True, ylim=[-250, 250], fill_feas=True)
        self.plot_f_leg(ax_RH, fx_RH, fy_RH, fz_RH, mfz_RH, N1, N2, N3, dt1, dt2, dt3, 'Right-hip leg', True, ylim=[-250, 250], fill_feas=True)

        plt.legend(edgecolor='black', loc='upper center', bbox_to_anchor=(-0.2, 3.0), ncol=4)

        plt.savefig(
            fig_name,
            bbox_inches="tight", 
            pad_inches=0.1
        )

        # self.plot_f_leg(ax_LF, fx_LF, fy_LF, fz_LF, mfz_LF, N1, N2, N3, dt1, dt2, dt3, 'Left-front leg', False, ylim=[-200, 200], fill_feas=True)
        # self.plot_f_leg(ax_LH, fx_LH, fy_LH, fz_LH, mfz_LH, N1, N2, N3, dt1, dt2, dt3, 'Left-hip leg', False, ylim=[-100, 100], fill_feas=True)
        # self.plot_f_leg(ax_RF, fx_RF, fy_RF, fz_RF, mfz_RF, N1, N2, N3, dt1, dt2, dt3, 'Right-front leg', True, ylim=[-200, 200], fill_feas=True)
        # self.plot_f_leg(ax_RH, fx_RH, fy_RH, fz_RH, mfz_RH, N1, N2, N3, dt1, dt2, dt3, 'Right-hip leg', True, ylim=[-100, 100], fill_feas=True)

        # ax_LF.legend_ = None
        # ax_LH.legend_ = None
        # ax_RF.legend_ = None
        # ax_RH.legend_ = None
        # plt.gca.legend_ = None

        # plt.savefig(
        #     'ylim_'+fig_name,
        #     bbox_inches="tight", 
        #     pad_inches=0.1
        # )
