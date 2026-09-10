"""
Excitability class of a conductance-based neuron: Hodgkin-Huxley vs Connor-Stevens.

Motivation
----------
The textbook Hodgkin-Huxley model has a *discontinuous* f-I curve: firing switches on
at a nonzero rate of roughly 50 Hz (Class 2, born at a subcritical Hopf bifurcation).
Cortical pyramidal neurons do not behave that way.  They fire regularly at a few Hz
and their f-I curve rises continuously from zero (Class 1, born at a SNIC).  The
tempting conclusion is that conductance-based modelling is too crude to reproduce real
cortical excitability.  That conclusion is wrong, and this experiment shows why: the
1952 parameter set describes one specific structure, the squid giant axon, and the
excitability class is set by the *balance of currents near threshold*, not by the
modelling framework.  Connor and Stevens changed that balance by adding one current,
the A-type potassium current I_A, whose inactivation gate acts as an amplifying
variable.  Does that single addition move the membrane from Class 2 to Class 1?

What it measures
----------------
For each model it sweeps the injected current I across rheobase, integrates the
membrane equations with forward Euler, discards the transient, and reads the steady
firing rate f from the interval spanned by the remaining spikes.  Rheobase is then
bracketed with a refined sweep, because the whole question lives in a narrow window of
current.  Two diagnostics matter: (i) whether f leaves zero continuously or jumps, and
(ii) the lowest sustained rate just above rheobase.  Expected result: HH jumps to about
50 Hz, Connor-Stevens comes up continuously and sustains a few Hz.

Usage
-----
    python fi_class1_class2.py                          # both sweeps + figure
    python fi_class1_class2.py --out fig.png            # choose the output file
    python fi_class1_class2.py --model hh               # print the HH f-I table
    python fi_class1_class2.py --model cs               # print the Connor-Stevens table
    python fi_class1_class2.py --npoints 120 --T 4000   # finer / longer sweep
"""

import argparse

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

Cm = 1.0                                  # uF/cm^2, both models


def _safe(num, den):
    """alpha_m and alpha_n have a removable 0/0 singularity; take the limit there."""
    small = np.abs(den) < 1e-9
    return np.where(small, 1.0, num / np.where(small, 1.0, den))


# ------------------------------------------------- model 1: Hodgkin-Huxley 1952 (squid)
HH = dict(gNa=120.0, gK=36.0, gL=0.3, gA=0.0,
          ENa=50.0, EK=-77.0, EL=-54.387, EA=-77.0, V0=-65.0,
          name="hh")


def hh_rates(V):
    am = _safe(0.1 * (V + 40), 1 - np.exp(-(V + 40) / 10))
    bm = 4.0 * np.exp(-(V + 65) / 18)
    ah = 0.07 * np.exp(-(V + 65) / 20)
    bh = 1.0 / (1 + np.exp(-(V + 35) / 10))
    an = _safe(0.01 * (V + 55), 1 - np.exp(-(V + 55) / 10))
    bn = 0.125 * np.exp(-(V + 65) / 80)
    return am, bm, ah, bh, an, bn


# ------------------------------------ model 2: Connor-Stevens (HH kinetics + A-current)
CS = dict(gNa=120.0, gK=20.0, gL=0.3, gA=47.7,
          ENa=55.0, EK=-72.0, EL=-17.0, EA=-75.0, V0=-68.0,
          name="cs")


def cs_rates(V):
    am = _safe(0.38 * (V + 29.7), 1 - np.exp(-(V + 29.7) / 10))
    bm = 15.2 * np.exp(-0.0556 * (V + 54.7))
    ah = 0.266 * np.exp(-0.05 * (V + 48.0))
    bh = 3.8 / (1 + np.exp(-0.1 * (V + 18.0)))
    an = _safe(0.02 * (V + 45.7), 1 - np.exp(-(V + 45.7) / 10))
    bn = 0.25 * np.exp(-0.0125 * (V + 55.7))
    return am, bm, ah, bh, an, bn


def a_inf(V):
    return (0.0761 * np.exp(0.0314 * (V + 94.22)) / (1 + np.exp(0.0346 * (V + 1.17)))) ** (1 / 3)


def tau_a(V):
    return 0.3632 + 1.158 / (1 + np.exp(0.0497 * (V + 55.96)))


def b_inf(V):
    return (1.0 / (1 + np.exp(0.0688 * (V + 53.3)))) ** 4


def tau_b(V):
    return 1.24 + 2.678 / (1 + np.exp(0.0624 * (V + 50.0)))


def steady(a, b):
    return a / (a + b)


def firing_rate(I, model=HH, T=2000.0, dt=0.005, skip=600.0):
    """Steady firing rate (Hz) under constant current I (uA/cm^2), by Euler integration.

    Vectorised over I: pass an array and get an array of rates back.
    """
    rates = cs_rates if model["name"] == "cs" else hh_rates
    gNa, gK, gL, gA = model["gNa"], model["gK"], model["gL"], model["gA"]
    ENa, EK, EL, EA = model["ENa"], model["EK"], model["EL"], model["EA"]

    I = np.atleast_1d(np.asarray(I, float))
    V = np.full_like(I, model["V0"])
    am, bm, ah, bh, an, bn = rates(V)
    m, h, n = steady(am, bm), steady(ah, bh), steady(an, bn)
    a, b = a_inf(V), b_inf(V)

    steps = int(T / dt)
    i_skip = int(skip / dt)
    count = np.zeros_like(I)
    t_first = np.full_like(I, np.nan)
    t_last = np.full_like(I, np.nan)

    for i in range(steps):
        t = i * dt
        v = V
        am, bm, ah, bh, an, bn = rates(v)
        INa = gNa * m ** 3 * h * (v - ENa)
        IK = gK * n ** 4 * (v - EK)
        IL = gL * (v - EL)
        IA = gA * a ** 3 * b * (v - EA) if gA else 0.0
        V = v + dt * (I - INa - IK - IL - IA) / Cm
        m = m + dt * (am * (1 - m) - bm * m)
        h = h + dt * (ah * (1 - h) - bh * h)
        n = n + dt * (an * (1 - n) - bn * n)
        if gA:
            a = a + dt * (a_inf(v) - a) / tau_a(v)
            b = b + dt * (b_inf(v) - b) / tau_b(v)

        if i >= i_skip:
            crossed = (v < 0.0) & (V >= 0.0)          # upward crossing of 0 mV
            if crossed.any():
                t_first = np.where(crossed & np.isnan(t_first), t, t_first)
                t_last = np.where(crossed, t, t_last)
                count = count + crossed

    span = t_last - t_first
    with np.errstate(invalid="ignore", divide="ignore"):
        f = np.where(count >= 2, 1000.0 * (count - 1) / span, 0.0)
    return np.nan_to_num(f)


def sweep(model, I_lo, I_hi, npoints, T):
    """Coarse sweep, then a refined sweep around the detected rheobase."""
    I = np.linspace(I_lo, I_hi, npoints)
    f = firing_rate(I, model=model, T=T)
    if not f.any():
        return I, f
    k = int(np.argmax(f > 0))
    lo = I[k - 1] if k else I_lo
    fine = np.linspace(lo, I[k], 25)
    I = np.concatenate([I[:k], fine, I[k + 1:]])
    f = np.concatenate([f[:k], firing_rate(fine, model=model, T=T), f[k + 1:]])
    return I, f


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--out", default="hh_connor_stevens_fi.png")
    p.add_argument("--model", choices=["hh", "cs"], default=None,
                   help="print one model's f-I table instead of drawing the figure")
    p.add_argument("--npoints", type=int, default=60)
    p.add_argument("--T", type=float, default=2000.0)
    args = p.parse_args()

    if args.model:
        model = CS if args.model == "cs" else HH
        I, f = sweep(model, 2.0, 25.0, args.npoints, args.T)
        for ii, ff in zip(I, f):
            print(f"I = {ii:7.4f}   f = {ff:7.2f} Hz")
        return

    I_hh, f_hh = sweep(HH, 2.0, 25.0, args.npoints, args.T)
    I_cs, f_cs = sweep(CS, 2.0, 25.0, args.npoints, args.T)

    fig, ax = plt.subplots(1, 2, figsize=(11, 4))
    for a_, I_, f_, title, sub in (
        (ax[0], I_hh, f_hh, "Hodgkin-Huxley 1952  (squid giant axon)",
         "Class 2: firing switches on at a nonzero rate"),
        (ax[1], I_cs, f_cs, r"Connor-Stevens  (same framework $+\ I_A$)",
         "Class 1: rate rises continuously from zero"),
    ):
        a_.plot(I_, f_, "-", lw=2, color="tab:blue")
        on = f_ > 0
        a_.plot(I_[on][0], f_[on][0], "o", ms=6, color="tab:red", zorder=5)
        a_.annotate(f"{f_[on][0]:.0f} Hz", (I_[on][0], f_[on][0]),
                    textcoords="offset points", xytext=(8, -4), color="tab:red")
        a_.set_title(f"{title}\n{sub}", fontsize=10)
        a_.set_xlabel(r"$I$  ($\mu$A/cm$^2$)")
        a_.set_ylabel(r"firing rate  $f$  (Hz)")
        a_.spines[["top", "right"]].set_visible(False)
        a_.set_ylim(bottom=0)
    fig.tight_layout()
    fig.savefig(args.out, dpi=150)

    for label, I_, f_ in (("HH            ", I_hh, f_hh), ("Connor-Stevens", I_cs, f_cs)):
        on = f_ > 0
        print(f"{label}: rheobase ~ {I_[on][0]:.3f} uA/cm^2, "
              f"lowest sustained rate {f_[on].min():.1f} Hz")
    print(f"wrote {args.out}")


if __name__ == "__main__":
    main()
