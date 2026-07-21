"""Figures for ch-biophys-03 (ion channels & gating, quantitative)."""
import numpy as np
import matplotlib.pyplot as plt
from style import *

OUT = "/tmp/build/ch-biophysics/"
rng = np.random.default_rng(7)

# ---- generic voltage-gated activation/inactivation (Boltzmann + bell tau) ----
def x_inf(V, Vh, k):          # k>0 activation, k<0 inactivation
    return 1.0/(1.0+np.exp(-(V-Vh)/k))

def tau_x(V, Vmax, tmin, tamp, sig):
    return tmin + tamp*np.exp(-((V-Vmax)/sig)**2)

# activation gate m: Vh=-25, k=+9 ; inactivation gate h: Vh=-45, k=-7
mVh, mk = -25.0, 9.0
hVh, hk = -45.0, -7.0

# ---------------------------------------------------------------
# Fig 1: steady-state activation/inactivation + time constants
# ---------------------------------------------------------------
def fig_gating_curves():
    V = np.linspace(-90, 40, 400)
    fig, axes = plt.subplots(1, 2, figsize=(9.6, 4.2))
    ax = axes[0]
    ax.plot(V, x_inf(V, mVh, mk), color=BLUE, label=r"activation  $m_\infty(V)$")
    ax.plot(V, x_inf(V, hVh, hk), color=RED, label=r"inactivation  $h_\infty(V)$")
    ax.axvline(mVh, color=BLUE, ls=":", lw=1); ax.axvline(hVh, color=RED, ls=":", lw=1)
    ax.fill_between(V, x_inf(V,mVh,mk)*x_inf(V,hVh,hk), color=PURPLE, alpha=0.18)
    ax.plot(V, x_inf(V,mVh,mk)*x_inf(V,hVh,hk), color=PURPLE, lw=1.6, ls="--",
            label=r"window  $m_\infty h_\infty$")
    ax.set_xlabel("membrane potential  $V$  (mV)")
    ax.set_ylabel("open probability")
    ax.set_title("(a) steady-state gating curves")
    ax.legend(loc="center left", fontsize=9.5)
    ax.set_ylim(-0.03, 1.05)

    ax = axes[1]
    ax.plot(V, tau_x(V, -20, 0.3, 3.2, 24), color=BLUE, label=r"$\tau_m(V)$")
    ax.plot(V, tau_x(V, -50, 1.5, 7.5, 20), color=RED, label=r"$\tau_h(V)$")
    ax.set_xlabel("membrane potential  $V$  (mV)")
    ax.set_ylabel("time constant  (ms)")
    ax.set_title("(b) voltage-dependent time constants")
    ax.legend(loc="upper right")
    fig.tight_layout()
    save(fig, OUT + "bio_gating_curves.png")

# ---------------------------------------------------------------
# Fig 2: simulated voltage clamp of an activating conductance g = gbar*m^p
# ---------------------------------------------------------------
def fig_voltage_clamp():
    gbar, E, p = 10.0, -80.0, 4     # delayed-rectifier-like K conductance
    dt, Tend = 0.02, 20.0
    steps = int(Tend/dt); t = np.arange(steps)*dt
    Vhold = -75.0
    steps_V = [-40, -20, 0, 20]
    fig, axes = plt.subplots(2, 1, figsize=(7.6, 6.4), sharex=True)
    cmap = [GREEN, BLUE, PURPLE, RED]
    for Vc, c in zip(steps_V, cmap):
        m = x_inf(Vhold, mVh, mk)*np.ones(steps)
        g = np.zeros(steps); I = np.zeros(steps)
        for k in range(steps):
            V = Vc if t[k] >= 2.0 else Vhold
            minf = x_inf(V, mVh, mk); tau = tau_x(V, -20, 0.3, 3.2, 24)
            if k < steps-1:
                m[k+1] = m[k] + dt*(minf-m[k])/tau
            g[k] = gbar*m[k]**p
            I[k] = g[k]*(V-E)
        axes[0].plot(t, g, color=c, label=fr"$V\!\to${Vc} mV")
        axes[1].plot(t, I, color=c)
    axes[0].axvline(2.0, color=LIGHT, ls=":", lw=1)
    axes[1].axvline(2.0, color=LIGHT, ls=":", lw=1)
    axes[0].set_ylabel(r"conductance  $g=\bar g\,m^4$")
    axes[0].set_title("Simulated voltage clamp of an activating conductance")
    axes[0].legend(loc="center right", fontsize=9.5)
    axes[1].set_ylabel("current  $I$")
    axes[1].set_xlabel("time (ms)")
    axes[1].text(2.2, axes[1].get_ylim()[1]*0.8, "voltage step", color=GRAY, fontsize=9)
    fig.tight_layout()
    save(fig, OUT + "bio_voltage_clamp.png")

# ---------------------------------------------------------------
# Fig 3: stochastic two-state channels -> ensemble -> deterministic limit
# ---------------------------------------------------------------
def fig_stochastic():
    # two-state C<->O, voltage-step turns activation on at t=2 ms
    dt, Tend = 0.05, 30.0
    steps = int(Tend/dt); t = np.arange(steps)*dt
    def rates(tt):
        # before step: mostly closed; after: alpha up
        if tt < 2.0:
            return 0.02, 1.0     # alpha (C->O), beta (O->C)
        return 1.0, 0.25
    minf = np.array([a/(a+b) for a,b in (rates(tt) for tt in t)])
    tau  = np.array([1.0/(a+b) for a,b in (rates(tt) for tt in t)])
    # deterministic
    m = np.zeros(steps); m[0] = 0.02/(0.02+1.0)
    for k in range(steps-1):
        a,b = rates(t[k]); m[k+1] = m[k] + dt*(a*(1-m[k]) - b*m[k])

    def simulate(N):
        state = (rng.random(N) < m[0]).astype(int)   # 0 closed 1 open
        frac = np.zeros(steps)
        traj0 = np.zeros(steps)
        for k in range(steps):
            frac[k] = state.mean()
            traj0[k] = state[0]
            a,b = rates(t[k])
            po_c = a*dt   # C->O prob
            po_o = b*dt   # O->C prob
            r = rng.random(N)
            newstate = state.copy()
            newstate[(state==0) & (r<po_c)] = 1
            newstate[(state==1) & (r<po_o)] = 0
            state = newstate
        return frac, traj0

    fig, axes = plt.subplots(1, 3, figsize=(11.4, 3.9))
    # panel 1: single channel traces
    ax = axes[0]
    for j, c in enumerate([BLUE, PURPLE, GREEN]):
        _, tr = simulate(1)
        ax.step(t, tr*0.8 + j*1.15, where="post", color=c, lw=1.3)
    ax.set_yticks([0.4, 1.55, 2.7]); ax.set_yticklabels(["ch 1","ch 2","ch 3"])
    ax.axvline(2.0, color=LIGHT, ls=":", lw=1)
    ax.set_xlabel("time (ms)")
    ax.set_title("(a) single channels: open/closed")
    ax.grid(axis="y", alpha=0)

    # panel 2: ensemble averages for increasing N
    ax = axes[1]
    for N, c in [(10, GREEN), (100, PURPLE), (1000, RED)]:
        frac, _ = simulate(N)
        ax.plot(t, frac, color=c, lw=1.4, alpha=0.9, label=f"N={N}")
    ax.plot(t, m, color="k", lw=2.4, ls="--", label=r"deterministic  $m(t)$")
    ax.axvline(2.0, color=LIGHT, ls=":", lw=1)
    ax.set_xlabel("time (ms)"); ax.set_ylabel("fraction open")
    ax.set_title("(b) ensemble average $\\to$ deterministic")
    ax.legend(loc="lower right", fontsize=9)

    # panel 3: noise amplitude ~ 1/sqrt(N)
    ax = axes[2]
    Ns = np.array([5, 10, 20, 50, 100, 300, 1000, 3000])
    sd = []
    for N in Ns:
        frac, _ = simulate(int(N))
        sd.append(np.std(frac[t>15]))   # steady-state fluctuation
    sd = np.array(sd)
    ax.loglog(Ns, sd, "o", color=BLUE, ms=7)
    ref = sd[0]*np.sqrt(Ns[0])/np.sqrt(Ns)
    ax.loglog(Ns, ref, color=RED, ls="--", label=r"$\propto 1/\sqrt{N}$")
    ax.set_xlabel("number of channels  $N$")
    ax.set_ylabel("std. of open fraction")
    ax.set_title("(c) relative noise shrinks as $1/\\sqrt{N}$")
    ax.legend(loc="upper right")
    fig.tight_layout()
    save(fig, OUT + "bio_stochastic_channels.png")

if __name__ == "__main__":
    fig_gating_curves()
    fig_voltage_clamp()
    fig_stochastic()
