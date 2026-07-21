"""Figures for ch-biophys-02 (computational membrane: Nernst / GHK / RC)."""
import numpy as np
import matplotlib.pyplot as plt
from style import *

OUT = "/tmp/build/ch-biophysics/"

# ---- physical constants ----
R = 8.314          # J/(K mol)
F = 96485.0        # C/mol
T = 310.0          # K (body temp)
RT_F = R * T / F * 1e3   # mV  (~26.71 mV)

# concentrations (mM): (out, in, z)
ions = {
    "K+":  (4.0,   140.0,  +1),
    "Na+": (145.0, 12.0,   +1),
    "Cl-": (120.0, 4.0,    -1),
    "Ca2+":(1.8,   1e-4,   +2),
}

def nernst(o, i, z):
    return RT_F / z * np.log(o / i)

# ---------------------------------------------------------------
# Fig 1: Nernst equilibrium potentials (bar chart)
# ---------------------------------------------------------------
def fig_nernst():
    names = list(ions.keys())
    labels = [r"$\mathrm{K}^+$", r"$\mathrm{Na}^+$", r"$\mathrm{Cl}^-$", r"$\mathrm{Ca}^{2+}$"]
    E = [nernst(*ions[n]) for n in names]
    colors = [BLUE if e < 0 else RED for e in E]
    fig, ax = plt.subplots(figsize=(7.3, 4.3))
    bars = ax.bar(labels, E, color=colors, width=0.6, edgecolor="#333", linewidth=0.8)
    ax.axhline(-65, color=GRAY, ls="--", lw=1.5)
    ax.text(3.35, -61, "resting  $V_m \\approx -65$", color=GRAY, ha="right", fontsize=10)
    ax.axhline(0, color="#444", lw=1.0)
    for b, e in zip(bars, E):
        va = "bottom" if e >= 0 else "top"
        off = 3 if e >= 0 else -3
        ax.text(b.get_x()+b.get_width()/2, e+off, f"{e:.0f}", ha="center", va=va, fontsize=10.5)
    ax.set_ylabel("equilibrium potential  $E_X$  (mV)")
    ax.set_title("Nernst equilibrium potentials (37 °C)")
    ax.set_ylim(-115, 155)
    save(fig, OUT + "bio_nernst_bars.png")

# ---------------------------------------------------------------
# Fig 2: Goldman/GHK resting potential vs P_Na/P_K, and during a
#        transient Na-permeability increase (an action-potential-like swing)
# ---------------------------------------------------------------
def ghk_V(PK, PNa, PCl):
    Ko, Ki, _ = ions["K+"]; Nao, Nai, _ = ions["Na+"]; Clo, Cli, _ = ions["Cl-"]
    num = PK*Ko + PNa*Nao + PCl*Cli
    den = PK*Ki + PNa*Nai + PCl*Clo
    return RT_F * np.log(num / den)

def fig_goldman():
    fig, axes = plt.subplots(1, 2, figsize=(9.6, 4.2))

    # (a) V_rest vs permeability ratio
    r = np.logspace(-3, 1.3, 400)          # PNa/PK
    PCl_ratio = 0.45
    V = np.array([ghk_V(1.0, ri, PCl_ratio) for ri in r])
    ax = axes[0]
    ax.semilogx(r, V, color=BLUE, lw=2.3)
    EK = nernst(*ions["K+"]); ENa = nernst(*ions["Na+"])
    ax.axhline(EK, color=GRAY, ls=":", lw=1.4); ax.text(2e-3, EK+3, "$E_K$", color=GRAY, fontsize=10)
    ax.axhline(ENa, color=GRAY, ls=":", lw=1.4); ax.text(2e-3, ENa-9, "$E_{Na}$", color=GRAY, fontsize=10)
    ax.plot(0.04, ghk_V(1,0.04,PCl_ratio), "o", color=RED, ms=8, zorder=5)
    ax.annotate("rest\n$P_{Na}/P_K\\approx0.04$", (0.04, ghk_V(1,0.04,PCl_ratio)),
                xytext=(0.0016, -30), fontsize=9.5, color=RED,
                arrowprops=dict(arrowstyle="->", color=RED, lw=1.2))
    ax.set_xlabel("permeability ratio  $P_{Na}/P_{K}$")
    ax.set_ylabel("membrane potential  $V_m$  (mV)")
    ax.set_title("(a) Goldman potential vs. $P_{Na}/P_K$")
    ax.set_ylim(-105, 75)

    # (b) transient: PNa/PK jumps up briefly -> V swings toward ENa
    t = np.linspace(0, 12, 1200)
    base = 0.04
    peak = 20.0
    # a fast rise, slower fall bump in PNa/PK
    bump = peak*np.exp(-((t-3.0)/0.35)**2) * (t>=3.0) + peak*np.exp(-((t-3.0)/1.1)**2)*(t<3.0)
    # simpler: alpha-like transient
    bump = np.where(t>=2.5, peak*((t-2.5)/0.4)*np.exp(1-(t-2.5)/0.4), 0.0)
    ratio = base + bump
    Vt = np.array([ghk_V(1.0, rr, PCl_ratio) for rr in ratio])
    ax = axes[1]
    ax.plot(t, Vt, color=PURPLE, lw=2.3)
    ax.axhline(EK, color=GRAY, ls=":", lw=1.2); ax.axhline(ENa, color=GRAY, ls=":", lw=1.2)
    ax.text(11.7, EK+3, "$E_K$", color=GRAY, ha="right", fontsize=10)
    ax.text(11.7, ENa-9, "$E_{Na}$", color=GRAY, ha="right", fontsize=10)
    ax.set_xlabel("time (a.u.)")
    ax.set_ylabel("$V_m$  (mV)")
    ax.set_title("(b) transient rise of $P_{Na}$ swings $V_m$ up")
    ax.set_ylim(-105, 75)
    fig.tight_layout()
    save(fig, OUT + "bio_goldman_sweep.png")

# ---------------------------------------------------------------
# Fig 3: RC membrane step response (charging), several time constants
# ---------------------------------------------------------------
def fig_rc_step():
    EL = -65.0
    Rm = 100.0        # MΩ
    I = 0.30          # nA  -> steady displacement I*Rm = 30 mV
    fig, axes = plt.subplots(1, 2, figsize=(9.6, 4.2))

    # (a) charging for different tau (different Cm)
    ax = axes[0]
    tvec = np.linspace(0, 80, 2000)
    for tau, c in [(5.0, GREEN), (10.0, BLUE), (20.0, PURPLE)]:
        Cm = tau / Rm * 1e3   # pF ( tau[ms]=Rm[MΩ]*Cm[pF]*1e-3 )
        Vinf = EL + I*Rm
        V = EL + (Vinf-EL)*(1-np.exp(-tvec/tau))
        ax.plot(tvec, V, color=c, label=fr"$\tau_m={tau:.0f}$ ms")
        ax.plot(tau, EL+(Vinf-EL)*(1-np.exp(-1)), "o", color=c, ms=6)
    Vinf = EL + I*Rm
    ax.axhline(Vinf, color=GRAY, ls="--", lw=1.2)
    ax.text(78, Vinf+1.2, r"$V_\infty = E_L + I R_m$", color=GRAY, ha="right", fontsize=10)
    ax.axhline(EL, color=LIGHT, ls=":", lw=1.2)
    ax.text(2, EL-3.4, "rest $E_L$", color=GRAY, fontsize=9.5)
    ax.set_xlabel("time (ms)"); ax.set_ylabel("$V_m$  (mV)")
    ax.set_title("(a) charging to a current step")
    ax.legend(loc="lower right")

    # (b) numerical Euler vs exact (short integration demo)
    ax = axes[1]
    tau = 10.0; Cm = tau/Rm*1e3
    dt = 1.0; Tend = 80
    steps = int(Tend/dt)
    tv = np.arange(steps)*dt
    V = np.zeros(steps); V[0] = EL
    def Iext(tt): return I if tt>=10 else 0.0
    for k in range(steps-1):
        V[k+1] = V[k] + dt*((EL - V[k])/Rm + Iext(tv[k]))/Cm*1e3  # keep mV/ms consistent
    # exact solution
    tf = np.linspace(0, Tend, 1000)
    Vex = np.where(tf<10, EL, EL + I*Rm*(1-np.exp(-(tf-10)/tau)))
    ax.plot(tf, Vex, color=GRAY, lw=3, alpha=0.6, label="exact")
    ax.plot(tv, V, color=BLUE, marker="o", ms=3.5, lw=1.2, label=r"Euler $\Delta t=1$ ms")
    ax.axvline(10, color=LIGHT, ls=":", lw=1)
    ax.text(11, -63, "step on", fontsize=9, color=GRAY)
    ax.set_xlabel("time (ms)"); ax.set_ylabel("$V_m$  (mV)")
    ax.set_title("(b) integrating the RC equation from scratch")
    ax.legend(loc="lower right")
    fig.tight_layout()
    save(fig, OUT + "bio_rc_step.png")

# ---------------------------------------------------------------
# Fig 4: strength-duration curve (rheobase & chronaxie) for RC + threshold
# ---------------------------------------------------------------
def fig_strength_duration():
    EL = -65.0; Vth = -50.0; Rm = 100.0; tau = 10.0
    Vth_disp = Vth - EL              # 15 mV
    Irheo = Vth_disp / Rm            # nA  (0.15 nA)
    Tvec = np.linspace(0.5, 60, 500)
    Ith = Irheo / (1 - np.exp(-Tvec/tau))
    chronaxie = tau*np.log(2)
    fig, ax = plt.subplots(figsize=(7.3, 4.3))
    ax.plot(Tvec, Ith, color=BLUE, lw=2.4)
    ax.axhline(Irheo, color=RED, ls="--", lw=1.6)
    ax.text(58, Irheo+0.012, "rheobase $I_{rh}$", color=RED, ha="right", fontsize=10.5)
    ax.axhline(2*Irheo, color=GRAY, ls=":", lw=1.2)
    ax.plot(chronaxie, 2*Irheo, "o", color=PURPLE, ms=8, zorder=5)
    ax.annotate(f"chronaxie\n$\\tau\\ln 2 \\approx {chronaxie:.1f}$ ms", (chronaxie, 2*Irheo),
                xytext=(18, 0.42), fontsize=10, color=PURPLE,
                arrowprops=dict(arrowstyle="->", color=PURPLE, lw=1.2))
    ax.set_xlabel("pulse duration  $T$  (ms)")
    ax.set_ylabel("threshold current  $I_{th}$  (nA)")
    ax.set_title("Strength–duration curve of an RC membrane")
    ax.set_ylim(0, 0.75)
    save(fig, OUT + "bio_strength_duration.png")

if __name__ == "__main__":
    print("RT/F =", round(RT_F,2), "mV")
    for n in ions: print(n, "E =", round(nernst(*ions[n]),1), "mV")
    fig_nernst()
    fig_goldman()
    fig_rc_step()
    fig_strength_duration()
