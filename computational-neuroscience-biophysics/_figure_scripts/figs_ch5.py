"""Figures for ch-biophys-05 (synapses & synaptic transmission)."""
import numpy as np
import matplotlib.pyplot as plt
from style import *

OUT = "/tmp/build/ch-biophysics/"

# ---------------------------------------------------------------
# Fig 1: synaptic conductance waveforms
# ---------------------------------------------------------------
def fig_conductance():
    t = np.linspace(0, 30, 1500)
    # single exponential
    td = 5.0
    g_exp = np.exp(-t/td)
    # alpha function, peak at tau
    tau = 2.0
    g_alpha = (t/tau)*np.exp(1 - t/tau)
    # double exponential (rise tr, decay td2), normalized to peak 1
    tr, td2 = 0.5, 5.0
    raw = np.exp(-t/td2) - np.exp(-t/tr)
    g_dexp = raw/raw.max()
    fig, ax = plt.subplots(figsize=(7.6, 4.4))
    ax.plot(t, g_exp,  color=BLUE,   label=r"single exp:  $e^{-t/\tau_d}$")
    ax.plot(t, g_alpha,color=GREEN,  label=r"alpha:  $(t/\tau)\,e^{1-t/\tau}$")
    ax.plot(t, g_dexp, color=PURPLE, label=r"double exp:  $e^{-t/\tau_d}-e^{-t/\tau_r}$")
    ax.set_xlabel("time since release  (ms)")
    ax.set_ylabel(r"normalized conductance  $g_{syn}(t)/\bar g$")
    ax.set_title("Synaptic conductance waveforms")
    ax.legend(loc="upper right")
    ax.set_ylim(-0.03, 1.08)
    save(fig, OUT + "bio_syn_conductance.png")

# ---------------------------------------------------------------
# Fig 2: EPSP vs IPSP (and shunting) through an RC membrane
# ---------------------------------------------------------------
def fig_epsp_ipsp():
    EL, gL, C = -65.0, 0.1, 1.0     # tau = C/gL = 10 ms
    dt, Tend = 0.05, 120.0
    steps = int(Tend/dt); t = np.arange(steps)*dt
    tr, td = 0.5, 5.0
    def gsyn(tt, t0, amp):
        s = tt - t0
        if s < 0: return 0.0
        raw = np.exp(-s/td) - np.exp(-s/tr)
        return amp*raw/0.63      # ~normalize peak

    def run(Esyn, amp, t0=10.0):
        V = np.zeros(steps); V[0] = EL
        for k in range(steps-1):
            g = gsyn(t[k], t0, amp)
            V[k+1] = V[k] + dt*(-gL*(V[k]-EL) - g*(V[k]-Esyn))/C
        return V

    fig, axes = plt.subplots(1, 2, figsize=(9.8, 4.2))
    ax = axes[0]
    ax.plot(t, run(0.0, 0.05),  color=RED,   label="excitatory  $E_{syn}=0$ (EPSP)")
    ax.plot(t, run(-80.0, 0.05),color=BLUE,  label="inhibitory  $E_{syn}=-80$ (IPSP)")
    ax.axhline(EL, color=GRAY, ls=":", lw=1.2); ax.text(115, EL+0.4, "rest", color=GRAY, ha="right", fontsize=9.5)
    ax.set_xlabel("time (ms)"); ax.set_ylabel("$V_m$  (mV)")
    ax.set_title("(a) EPSP and IPSP")
    ax.legend(loc="upper right", fontsize=9.5); ax.set_xlim(0, 80)

    # (b) reversal potential: same synapse from different holding potentials
    ax = axes[1]
    for Vhold, c in [(-80, BLUE), (-65, GREEN), (-40, PURPLE), (-10, RED)]:
        EL2 = Vhold
        V = np.zeros(steps); V[0] = EL2
        for k in range(steps-1):
            g = gsyn(t[k], 10.0, 0.05)
            V[k+1] = V[k] + dt*(-gL*(V[k]-EL2) - g*(V[k]-0.0))/C
        ax.plot(t, V-EL2, color=c, label=f"hold {Vhold} mV")
    ax.axhline(0, color="#444", lw=1)
    ax.set_xlabel("time (ms)"); ax.set_ylabel(r"$V_m - V_{hold}$  (mV)")
    ax.set_title(r"(b) same synapse, different holding $V$")
    ax.legend(loc="upper right", fontsize=9); ax.set_xlim(0, 80)
    fig.tight_layout()
    save(fig, OUT + "bio_epsp_ipsp.png")

# ---------------------------------------------------------------
# Fig 3: NMDA Mg block
# ---------------------------------------------------------------
def fig_nmda():
    V = np.linspace(-90, 40, 400)
    Mg = 1.0
    B = 1.0/(1.0 + np.exp(-0.062*V)*Mg/3.57)
    fig, axes = plt.subplots(1, 2, figsize=(9.6, 4.2))
    ax = axes[0]
    for mg, c in [(0.0, GRAY),(0.5, GREEN),(1.0, BLUE),(2.0, PURPLE)]:
        Bm = 1.0/(1.0 + np.exp(-0.062*V)*mg/3.57)
        ax.plot(V, Bm, color=c, label=f"[Mg]={mg} mM")
    ax.set_xlabel("membrane potential  $V$  (mV)")
    ax.set_ylabel(r"fraction unblocked  $B(V)$")
    ax.set_title("(a) voltage-dependent Mg$^{2+}$ block")
    ax.legend(loc="upper left", fontsize=9)

    ax = axes[1]
    E = 0.0
    I_ampa = 1.0*(V-E)
    I_nmda = 1.0*B*(V-E)
    ax.plot(V, I_ampa, color=GREEN, label="AMPA  $g(V-E)$")
    ax.plot(V, I_nmda, color=RED,   label="NMDA  $g\\,B(V)(V-E)$")
    ax.axhline(0, color="#444", lw=1); ax.axvline(0, color=LIGHT, lw=1, ls=":")
    ax.set_xlabel("membrane potential  $V$  (mV)")
    ax.set_ylabel("synaptic current  $I$  (a.u.)")
    ax.set_title("(b) current–voltage relations")
    ax.legend(loc="upper left", fontsize=9.5)
    fig.tight_layout()
    save(fig, OUT + "bio_nmda_mgblock.png")

# ---------------------------------------------------------------
# Fig 4: short-term plasticity (Tsodyks-Markram)
# ---------------------------------------------------------------
def fig_stp():
    spikes = np.arange(8)*50.0 + 20.0   # 20 Hz train
    Tend = 480.0
    def tm_model(U, tau_rec, tau_fac):
        R, u = 1.0, U
        amps, times = [], []
        last = None
        for ts in spikes:
            if last is not None:
                dt = ts - last
                R = 1 - (1-R)*np.exp(-dt/tau_rec)
                if tau_fac > 0:
                    u = U + (u-U)*np.exp(-dt/tau_fac)
            if tau_fac > 0:
                u = u + U*(1-u)      # facilitation increment at spike
            A = u*R
            amps.append(A); times.append(ts)
            R = R - u*R              # deplete
            last = ts
        return np.array(times), np.array(amps)

    fig, ax = plt.subplots(figsize=(8.4, 4.4))
    for (U, trec, tfac), name, c in [
        ((0.6, 800.0, 0.0), "depressing (U=0.6, high release)", BLUE),
        ((0.15, 100.0, 500.0), "facilitating (U=0.15, low release)", RED)]:
        tm, am = tm_model(U, trec, tfac)
        am = am/am[0]
        ax.plot(tm, am, "o-", color=c, ms=7, label=name)
        for x, y in zip(tm, am):
            ax.plot([x, x], [0, y], color=c, lw=1, alpha=0.35)
    ax.set_xlabel("time (ms)")
    ax.set_ylabel("EPSC amplitude (norm. to first)")
    ax.set_title("Short-term plasticity for a 20 Hz spike train")
    ax.legend(loc="upper right", fontsize=9.5)
    ax.set_ylim(0, 2.6)
    save(fig, OUT + "bio_stp.png")

# ---------------------------------------------------------------
# Fig 5: temporal & spatial summation
# ---------------------------------------------------------------
def fig_summation():
    EL, gL, C, Vth = -65.0, 0.1, 1.0, -50.0
    dt, Tend = 0.05, 200.0
    steps = int(Tend/dt); t = np.arange(steps)*dt
    tr, td = 0.5, 6.0
    def gkernel(s):
        return np.where(s>=0, (np.exp(-s/td)-np.exp(-s/tr))/0.6, 0.0)
    def run(events, amp):
        # events: list of (t0) ; each excitatory synapse E=0
        V = np.zeros(steps); V[0] = EL
        for k in range(steps-1):
            g = sum(amp*gkernel(t[k]-t0) for t0 in events)
            V[k+1] = V[k] + dt*(-gL*(V[k]-EL) - g*(V[k]-0.0))/C
        return V
    fig, axes = plt.subplots(1, 2, figsize=(9.8, 4.2))
    # temporal: same synapse repeated at short vs long interval
    ax = axes[0]
    ax.plot(t, run([20], 0.04), color=GRAY, label="single input")
    ax.plot(t, run([20,35,50,65], 0.04), color=BLUE, label="fast train (sums)")
    ax.plot(t, run([20,70,120,170], 0.04), color=GREEN, label="slow train (no sum)")
    ax.axhline(Vth, color=RED, ls="--", lw=1.4); ax.text(198, Vth+0.5, "threshold", color=RED, ha="right", fontsize=9.5)
    ax.set_xlabel("time (ms)"); ax.set_ylabel("$V_m$  (mV)")
    ax.set_title("(a) temporal summation")
    ax.legend(loc="upper right", fontsize=9)
    # spatial: increasing number of simultaneous synapses
    ax = axes[1]
    for nsyn, c in [(1, GRAY),(3, GREEN),(6, BLUE),(10, PURPLE)]:
        V = run([20]*nsyn, 0.04)   # nsyn synapses at once == amp*nsyn
        ax.plot(t, V, color=c, label=f"{nsyn} synapses")
    ax.axhline(Vth, color=RED, ls="--", lw=1.4); ax.text(198, Vth+0.5, "threshold", color=RED, ha="right", fontsize=9.5)
    ax.set_xlabel("time (ms)"); ax.set_ylabel("$V_m$  (mV)")
    ax.set_title("(b) spatial summation")
    ax.legend(loc="upper right", fontsize=9); ax.set_xlim(0, 80)
    fig.tight_layout()
    save(fig, OUT + "bio_summation.png")

if __name__ == "__main__":
    fig_conductance()
    fig_epsp_ipsp()
    fig_nmda()
    fig_stp()
    fig_summation()
