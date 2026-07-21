"""Figures for ch-biophys-04 (cable theory & dendrites)."""
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, FancyArrow
from style import *

OUT = "/tmp/build/ch-biophysics/"

# ---------------------------------------------------------------
# Fig 1: RC-ladder schematic of a passive cable
# ---------------------------------------------------------------
def fig_ladder():
    fig, ax = plt.subplots(figsize=(9.2, 3.7))
    n = 5
    x0, dx = 0.6, 1.9
    ytop, ybot = 2.2, 0.0
    # extracellular rail
    ax.plot([x0-0.5, x0+dx*(n-1)+0.5], [ybot, ybot], color="#444", lw=2)
    # intracellular rail with axial resistors (boxes)
    for i in range(n):
        x = x0 + i*dx
        # node
        ax.plot(x, ytop, "o", color="#444", ms=5, zorder=5)
        # membrane branch: cm and rm in PARALLEL from node to ground
        ax.plot([x, x], [ytop, ytop-0.28], color="#444", lw=1.6)     # short stem
        xc, xr = x-0.32, x+0.32                                       # two parallel legs
        ax.plot([xc, xr], [ytop-0.28, ytop-0.28], color="#444", lw=1.6)  # top tie
        ax.plot([xc, xr], [ybot+0.28, ybot+0.28], color="#444", lw=1.6)  # bottom tie
        ax.plot([x, x], [ybot+0.28, ybot], color="#444", lw=1.6)     # short stem to rail
        # capacitor leg (left)
        ax.plot([xc, xc], [ytop-0.28, ytop-0.62], color="#444", lw=1.6)
        ax.plot([xc-0.16, xc+0.16], [ytop-0.62, ytop-0.62], color=BLUE, lw=2.4)
        ax.plot([xc-0.16, xc+0.16], [ytop-0.78, ytop-0.78], color=BLUE, lw=2.4)
        ax.plot([xc, xc], [ytop-0.78, ybot+0.28], color="#444", lw=1.6)
        ax.text(xc-0.5, ytop-0.75, r"$c_m$", color=BLUE, fontsize=10)
        # resistor leg (right)
        ax.plot([xr, xr], [ytop-0.28, ytop-0.5], color="#444", lw=1.6)
        ax.add_patch(Rectangle((xr-0.12, ytop-1.05), 0.24, 0.55,
                     fill=True, facecolor="white", edgecolor=GREEN, lw=1.8))
        ax.plot([xr, xr], [ytop-1.05, ybot+0.28], color="#444", lw=1.6)
        ax.text(xr+0.2, ytop-0.82, r"$r_m$", color=GREEN, fontsize=10)
    # axial resistors between nodes
    for i in range(n-1):
        x = x0 + i*dx
        xr = x + dx*0.5
        ax.add_patch(Rectangle((xr-0.28, ytop-0.13), 0.56, 0.26,
                     fill=True, facecolor="white", edgecolor=RED, lw=1.8))
        ax.text(xr, ytop+0.28, r"$r_a$", color=RED, ha="center", fontsize=11)
        ax.plot([x+0.06, xr-0.28], [ytop, ytop], color="#444", lw=1.6)
        ax.plot([xr+0.28, x+dx-0.06], [ytop, ytop], color="#444", lw=1.6)
    # injected current arrow at node 0
    ax.annotate("", xy=(x0, ytop+0.05), xytext=(x0-0.9, ytop+0.05),
                arrowprops=dict(arrowstyle="-|>", color=PURPLE, lw=2.2))
    ax.text(x0-1.0, ytop+0.18, r"$I$", color=PURPLE, fontsize=13)
    ax.text(x0+dx*(n-1)+0.7, ytop, r"$\cdots$", fontsize=16, va="center")
    # x-axis label
    ax.annotate("", xy=(x0+dx*(n-1)+0.2, -0.75), xytext=(x0-0.3, -0.75),
                arrowprops=dict(arrowstyle="->", color=GRAY, lw=1.4))
    ax.text(x0+dx*(n-1)/2, -1.05, "distance along the cable  $x$", color=GRAY, ha="center", fontsize=11)
    ax.text(x0-0.6, ybot-0.35, "extracellular (ground)", color=GRAY, fontsize=9)
    ax.set_xlim(x0-1.3, x0+dx*(n-1)+1.3)
    ax.set_ylim(-1.3, 2.9)
    ax.axis("off")
    ax.set_title("A passive cable as a ladder of RC compartments")
    save(fig, OUT + "bio_cable_ladder.png")

# ---------------------------------------------------------------
# Fig 2: steady-state decay along the cable
# ---------------------------------------------------------------
def fig_steadystate():
    x = np.linspace(0, 5, 400)   # mm
    fig, axes = plt.subplots(1, 2, figsize=(9.6, 4.2))
    ax = axes[0]
    for lam, c in [(0.5, GREEN), (1.0, BLUE), (2.0, PURPLE)]:
        ax.plot(x, np.exp(-x/lam), color=c, label=fr"$\lambda={lam}$ mm")
        ax.plot(lam, np.exp(-1), "o", color=c, ms=6)
    ax.axhline(np.exp(-1), color=GRAY, ls=":", lw=1.2)
    ax.text(4.9, np.exp(-1)+0.02, r"$1/e\approx0.37$", color=GRAY, ha="right", fontsize=10)
    ax.set_xlabel("distance  $x$  (mm)")
    ax.set_ylabel(r"$V(x)/V_0$")
    ax.set_title(r"(a) semi-infinite cable:  $V=V_0 e^{-x/\lambda}$")
    ax.legend(loc="upper right")

    ax = axes[1]
    lam = 1.0
    for L, c in [(0.5, GREEN), (1.0, BLUE), (3.0, PURPLE)]:
        xx = np.linspace(0, L, 200)
        V = np.cosh((L-xx)/lam)/np.cosh(L/lam)
        ax.plot(xx, V, color=c, label=fr"$L={L}$ mm (sealed end)")
    ax.plot(x, np.exp(-x/lam), color=GRAY, ls="--", lw=1.5, label=r"infinite cable")
    ax.set_xlabel("distance  $x$  (mm)")
    ax.set_ylabel(r"$V(x)/V_0$")
    ax.set_title(r"(b) finite cable with sealed end ($\lambda=1$)")
    ax.legend(loc="upper right", fontsize=9)
    ax.set_xlim(0, 3.2)
    fig.tight_layout()
    save(fig, OUT + "bio_cable_steadystate.png")

# ---------------------------------------------------------------
# compartmental simulation of the discretized cable equation
#   dV/dt = (1/tau)[ lambda^2 V_xx - V ] + I/Cm
# ---------------------------------------------------------------
def simulate_cable(N=60, Lmm=6.0, lam=1.0, tau=10.0, dt=0.02, Tend=60.0,
                   inject=None):
    dx = Lmm/(N-1)
    coeff = (lam**2)/(dx**2)
    steps = int(Tend/dt)
    V = np.zeros(N)
    rec = np.zeros((steps, N))
    for k in range(steps):
        Vxx = np.zeros(N)
        Vxx[1:-1] = V[2:] - 2*V[1:-1] + V[:-2]
        Vxx[0]  = 2*(V[1]-V[0])     # sealed end (mirror)
        Vxx[-1] = 2*(V[-2]-V[-1])
        I = inject(k*dt) if inject else np.zeros(N)
        V = V + dt*((coeff*Vxx - V)/tau) + dt*I
        rec[k] = V
    return rec, dx, np.arange(steps)*dt

def fig_compartment():
    N = 60; Lmm = 6.0; lam = 1.0
    # steady current injected at x=0
    Iamp = np.zeros(N); Iamp[0] = 1.0
    rec, dx, t = simulate_cable(N=N, Lmm=Lmm, lam=lam, inject=lambda tt: Iamp)
    x = np.arange(N)*dx
    fig, axes = plt.subplots(1, 2, figsize=(9.8, 4.2))
    # (a) space-time heatmap
    ax = axes[0]
    im = ax.pcolormesh(x, t, rec, shading="auto", cmap="viridis")
    ax.set_xlabel("distance  $x$  (mm)"); ax.set_ylabel("time (ms)")
    ax.set_title("(a) voltage spread in space & time")
    cb = fig.colorbar(im, ax=ax); cb.set_label("$V$ (a.u.)")
    # (b) steady-state profile vs theory
    ax = axes[1]
    Vss = rec[-1]
    ax.plot(x, Vss/Vss[0], "o", color=BLUE, ms=4, label="simulation (steady state)")
    ax.plot(x, np.cosh((Lmm-x)/lam)/np.cosh(Lmm/lam), color=RED, ls="--",
            label=r"theory $\cosh$ (sealed end)")
    ax.set_xlabel("distance  $x$  (mm)"); ax.set_ylabel(r"$V(x)/V(0)$")
    ax.set_title(r"(b) steady profile matches cable theory")
    ax.legend(loc="upper right", fontsize=9.5)
    fig.tight_layout()
    save(fig, OUT + "bio_cable_compartment.png")

def fig_attenuation():
    # transient synaptic current injected at a distal site; record dendrite vs soma
    N = 60; Lmm = 6.0; lam = 1.0; tau = 10.0
    dx = Lmm/(N-1)
    site = 40   # distal compartment
    soma = 0
    def alpha_I(tt, t0=5.0, ts=2.0, amp=8.0):
        I = np.zeros(N)
        if tt >= t0:
            s = (tt-t0)/ts
            I[site] = amp*s*np.exp(1-s)
        return I
    rec, dx, t = simulate_cable(N=N, Lmm=Lmm, lam=lam, tau=tau, dt=0.02, Tend=60.0,
                                inject=alpha_I)
    fig, axes = plt.subplots(1, 2, figsize=(9.6, 4.2))
    ax = axes[0]
    ax.plot(t, rec[:, site], color=PURPLE, label=f"input site  $x={site*dx:.1f}$ mm")
    ax.plot(t, rec[:, soma], color=BLUE, label="soma  $x=0$")
    ax.set_xlabel("time (ms)"); ax.set_ylabel("$V$ (a.u.)")
    ax.set_title("(a) EPSP is attenuated, delayed & slowed")
    ax.legend(loc="upper right"); ax.set_xlim(0, 45)
    # (b) peak amplitude vs distance of recording (input fixed at 'site')
    ax = axes[1]
    peak = rec.max(axis=0)
    x = np.arange(N)*dx
    ax.plot(x, peak/peak[site], color=BLUE, lw=2.2)
    ax.axvline(site*dx, color=PURPLE, ls=":", lw=1.5)
    ax.text(site*dx-0.15, 0.5, "input site", color=PURPLE, rotation=90, va="center", fontsize=9.5)
    ax.plot(0, peak[soma]/peak[site], "o", color=BLUE, ms=8)
    ax.text(0.15, peak[soma]/peak[site]+0.03, "soma", color=BLUE, fontsize=10)
    ax.set_xlabel("recording position  $x$  (mm)")
    ax.set_ylabel("peak EPSP (norm. to input site)")
    ax.set_title("(b) peak attenuation with distance")
    fig.tight_layout()
    save(fig, OUT + "bio_dendrite_attenuation.png")

if __name__ == "__main__":
    fig_ladder()
    fig_steadystate()
    fig_compartment()
    fig_attenuation()
