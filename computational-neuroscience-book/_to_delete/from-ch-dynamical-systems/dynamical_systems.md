---
title: "Dynamical Systems"
subtitle: "Fixed points, stability, and bifurcations — from one dimension to neural circuits"
---

# Dynamical Systems

> *"All models are wrong, but some models are useful."* — George Box

Most of computational neuroscience is, at bottom, the study of **differential equations**: a
membrane potential that relaxes back to rest, a pair of genes that repress each other, a
population of neurons that synchronises into an oscillation. In almost every case we cannot
write the solution down as a formula. The good news is that we usually do not need to. With a
handful of ideas — *fixed points*, *stability*, *phase portraits*, and *bifurcations* — we can
predict the qualitative behaviour of a model across whole ranges of parameters and initial
conditions, **without ever solving it in closed form**. This branch of mathematics is called
**dynamical systems theory**, and it is one of the most powerful and beautiful tools available
to a modeller.

This chapter builds the theory from the ground up. We start in one dimension, where everything
can be read off a single picture (the *phase line*). We then move to systems of two or more
equations, where the central new object is the *Jacobian matrix*. Finally we put the whole
machine to work on two classic neuroscience-flavoured models: a **bistable genetic switch**
(the basic mechanism behind cellular memory) and the **FitzHugh–Nagumo neuron** (the simplest
model that fires action potentials).

:::{admonition} What you will be able to do by the end
:class: tip
- Find the **fixed points** of a model and decide whether each one is **stable** or **unstable**.
- Draw and interpret a **phase line** (1D) and a **phase plane** (2D).
- Classify a two-dimensional fixed point as a **node, saddle, or spiral** using the **trace** and **determinant** of the Jacobian.
- Recognise the two bifurcations that matter most in neuroscience — the **saddle-node** and the **Hopf** — and explain what each one does to the dynamics.
- Read code that turns all of the above into reproducible figures (phase planes, bifurcation diagrams).
:::

---

## 1. Autonomous systems: the idea of "letting the equation tell you what happens"

A first-order ordinary differential equation (ODE) has the form

$$
\frac{dx}{dt} = f(x, t), \qquad x(t_0) = x_0 .
$$

When the right-hand side does **not** depend explicitly on time — that is, $f(x,t)=f(x)$ — we
call the equation **autonomous**:

$$
\frac{dx}{dt} = f(x), \qquad x(t_0) = x_0 .
$$

Autonomy is the single property that makes the dynamical-systems viewpoint possible. The slope
$dx/dt$ at any given state $x$ is *always the same number*, no matter when you arrive there. The
state space is therefore covered by a fixed "flow field": at every value of $x$ there is one
prescribed velocity, and the solution is just whatever curve follows those velocities.

Many neuroscience models are autonomous whenever their input is held constant. The
leaky-integrate-and-fire (LIF) neuron driven by a constant current, or the FitzHugh–Nagumo
model with $I_\text{ext}$ fixed, are both autonomous: dividing through by the membrane time
constant puts them in exactly the form above. (When the input varies in time the system becomes
*non-autonomous*; we return to that case in {ref}`sec-nonauto`.)

For a much fuller and more leisurely treatment than this chapter can give, Steven Strogatz's
*Nonlinear Dynamics and Chaos* {cite}`strogatz` is the standard and very readable reference, and
the appendices of Rosenbaum's *Modeling Neural Circuits Made Simple with Python*
{cite}`rosenbaum2022` are the source that this chapter most closely follows.

---

## 2. One dimension: fixed points, phase lines, and stability

### 2.1 Fixed points

A **fixed point** (also called an *equilibrium* or *steady state*) of $\dot x = f(x)$ is a value
$x^*$ at which the flow vanishes:

$$
f(x^*) = 0 .
$$

We write $x^*$ (with a star) for fixed points to keep them distinct from the initial condition
$x_0$. The defining property of a fixed point is wonderfully simple:

:::{important}
**If a solution starts at a fixed point, it stays there forever.**
:::

The reason is immediate. If $x(t_0) = x^*$ and $f(x^*) = 0$, then $\dot x(t_0) = f(x^*) = 0$: the
state has no velocity, so it cannot move; with zero velocity it stays put, so its velocity
remains zero, and so on. The constant function $x(t)=x^*$ satisfies the ODE exactly. Fixed
points are therefore the *skeleton* of the dynamics — the resting states, the memory states, the
"do-nothing" configurations of a model.

### 2.2 Solutions cannot cross fixed points, and the velocity never changes sign

Two facts make one-dimensional dynamics extremely predictable. Both rely only on $f$ being
continuous and on solutions being unique.

1. **A trajectory can never reach a fixed point in finite "interesting" motion.** Suppose
   $\dot x(t_0) = f(x_0) > 0$ (the state is moving right) and suppose it later arrived exactly at
   some fixed point $x_1$. From $x_1$ the only solution is the constant one — but a constant
   solution has zero velocity *everywhere*, contradicting $\dot x(t_0)>0$. So a rightward-moving
   solution can approach a fixed point but never actually attain it.

2. **The sign of $\dot x$ is constant in time.** If the velocity were positive at one moment and
   negative at another, then by the intermediate value theorem it would have to pass through zero
   at some intermediate time — i.e. the trajectory would momentarily sit at a fixed point, which
   (by point 1) freezes it forever and again contradicts the assumption.

Putting these together:

:::{important}
For an autonomous 1-D ODE, **the sign of $\dot x(t)$ never changes**. Each solution moves
monotonically, either always increasing or always decreasing, and is trapped between adjacent
fixed points.
:::

### 2.3 The phase line

These observations mean we can summarise *every* solution in a single one-dimensional picture
called the **phase line**. The recipe is:

1. Plot $f(x)$ and mark the fixed points where it crosses zero.
2. On the $x$-axis, draw a **right-pointing arrow** wherever $f(x) > 0$ (the state increases) and
   a **left-pointing arrow** wherever $f(x) < 0$ (the state decreases).

The arrows tell you where any initial condition will flow. {numref}`fig-phase-line` does this for
the textbook example

$$
\dot x = x^2 - 1 ,
$$

which has fixed points at $x^* = -1$ and $x^* = +1$.

```{figure} figures/phase_line.png
:name: fig-phase-line
:width: 80%

Phase line for $\dot x = x^2 - 1$. The curve $f(x)=x^2-1$ is shown in teal; open circles on the
lower axis are the two fixed points. To the left of $-1$ and to the right of $+1$, $f>0$ and the
flow points right; in between, $f<0$ and the flow points left. The arrows therefore converge on
$x^*=-1$ (stable, blue) and diverge from $x^*=+1$ (unstable, red).
```

Notice how the **whole story** is visible at a glance: any state starting below $-1$ is pushed up
toward $-1$; any state between $-1$ and $+1$ slides down toward $-1$; any state above $+1$ runs
off to infinity. We have completely characterised the long-term behaviour without solving the
equation.

```python
import numpy as np
import matplotlib.pyplot as plt

def phase_line(f, xrange=(-2.2, 2.2), n=400):
    """Plot f(x) together with directional arrows showing the 1-D flow."""
    x = np.linspace(*xrange, n)
    y = f(x)
    fig, ax = plt.subplots(figsize=(6, 3.2))
    ax.axhline(0, color="0.6", lw=1)
    ax.plot(x, y, lw=2)
    # mark sign-change points (fixed points) and draw flow arrows
    sign = np.sign(y)
    for i in np.where(np.diff(sign) != 0)[0]:
        ax.plot(0.5 * (x[i] + x[i + 1]), 0, "o", ms=9, mfc="white", mec="k")
    for a, b in zip(x[:-1:40], x[1::40]):
        s = np.sign(f(0.5 * (a + b)))
        ax.annotate("", xy=(0.5*(a+b) + 0.1*s, 0), xytext=(0.5*(a+b) - 0.1*s, 0),
                    arrowprops=dict(arrowstyle="-|>", color="0.3"))
    ax.set(xlabel="x", ylabel="f(x)")
    return ax

phase_line(lambda x: x**2 - 1)
plt.show()
```

### 2.4 Stability and the linearisation test

Look again at {numref}`fig-phase-line`. Around $x^*=-1$ both arrows point *inward*: a state that
starts nearby is pushed back. Around $x^*=+1$ both arrows point *outward*: the smallest
displacement is amplified. This is the difference between a **stable** and an **unstable** fixed
point.

:::{admonition} Definitions
:class: note
A fixed point $x^*$ is **(asymptotically) stable** if every solution that starts sufficiently
close to $x^*$ converges to it as $t\to\infty$. It is **unstable** if solutions can start
arbitrarily close and still fail to converge.
:::

There is a quick algebraic test. Near a stable point the flow $f$ changes from positive (push
right) to negative (push left) as $x$ increases through $x^*$ — so $f$ is *decreasing* there.
Near an unstable point the reverse happens. This is exactly the sign of the derivative
$f'(x^*)$:

:::{important}
**Linear stability test (1-D).** Let $x^*$ be a fixed point of $\dot x = f(x)$.
- If $f'(x^*) < 0$, then $x^*$ is **stable**.
- If $f'(x^*) > 0$, then $x^*$ is **unstable**.
:::

For our example $f'(x) = 2x$, so $f'(-1) = -2 < 0$ (stable) and $f'(+1) = +2 > 0$ (unstable),
confirming the picture. The intuition — *negative slope pulls you back, positive slope throws you
out* — is the one-dimensional seed of everything that follows; in higher dimensions the role of
$f'(x^*)$ is taken over by the eigenvalues of a matrix.

### 2.5 Borderline cases

The test is silent when $f'(x^*) = 0$, because the curve is tangent to the axis and the linear
term tells us nothing. Two flavours of borderline behaviour can occur, and in both the phase line
still settles the matter instantly.

- **Semi-stable.** For $f(x) = x^2$ the only fixed point is $x^*=0$, with $f'(0)=0$. Trajectories
  approach from the left ($f<0$) but flee on the right ($f>0$). It attracts on one side and
  repels on the other.
- **Neutrally stable.** For the trivial flow $f(x)=0$ *every* point is a fixed point. Nearby
  states neither approach nor leave — they simply stay where they are.

Note that $f'(x^*)=0$ does **not** automatically mean "semi-stable": for $f(x)=x^3$ the point
$x^*=0$ is genuinely unstable. The lesson is that **drawing the phase line is more reliable than
memorising rules.** When $f'(x^*)\neq 0$, the linear test is decisive; otherwise, look at the
sign of $f$ on each side.

### 2.6 Our first bifurcation: the saddle-node

The real payoff of this viewpoint is that we can ask how the dynamics *change* as a parameter is
varied. A parameter value at which the **number or stability of the fixed points changes** is
called a **bifurcation**. In one dimension, the most common is the **saddle-node** (also called
*fold* or *tangent*) bifurcation, in which two fixed points collide and annihilate. The minimal
example is

$$
\dot x = x^2 + a ,
$$

where $a$ is a control parameter ({numref}`fig-saddle-node`).

```{figure} figures/saddle_node.png
:name: fig-saddle-node

Saddle-node bifurcation of $\dot x = x^2 + a$. **For $a<0$** there are two fixed points,
$x^* = \pm\sqrt{-a}$: the lower one is stable (blue), the upper unstable (red). **At $a=0$** they
merge into a single half-stable point. **For $a>0$** the parabola lifts off the axis and there are
*no* fixed points at all — the flow is everywhere positive. The right-hand panel is the
**bifurcation diagram**: fixed-point location versus $a$, with the stable branch solid and the
unstable branch dashed. The two branches meet and vanish at $a=0$.
```

In neuroscience this is the canonical way a neuron switches *on*: as the input current crosses a
threshold, a stable "rest" state and a nearby unstable "threshold" state collide and disappear,
leaving the system no choice but to fire. We will see exactly this happen in the FitzHugh–Nagumo
model in {ref}`sec-fhn`.

:::{admonition} Exercises (1-D)
:class: seealso
1. **Solve-then-check.** For $\dot x = x^2 - 1$ with $x(0)=1$, what is the exact solution? (Hint:
   $x_0$ is a fixed point.) Confirm numerically with forward Euler, $dt=0.01$, $T=5$.
2. **Phase line.** Sketch the phase line of $\dot x = -(x+1)(x-1)(x-2)$. Find all three fixed
   points and classify each with the sign of $f'(x^*)$. Verify by simulating from several initial
   conditions.
3. **Cubic.** Show that $x^*=0$ is unstable for $\dot x = x^3$ even though $f'(0)=0$. Why does the
   linear test fail here, and why does the phase line still work?
4. **Find the bifurcation.** For $\dot x = x^2 + a$, classify the fixed points for $a=-1$, $a=0$,
   and $a=1$, and confirm that the bifurcation is at $a=0$.
:::

---

## 3. Higher dimensions: systems of ODEs

Almost every interesting model needs at least two variables — a voltage and a recovery current,
two mutually-repressing genes, an excitatory and an inhibitory population. We therefore extend
everything above to **systems** of autonomous ODEs,

$$
\frac{d\mathbf{u}}{dt} = \mathbf{F}(\mathbf{u}), \qquad
\mathbf{u}\in\mathbb{R}^n,\;\; \mathbf{F}:\mathbb{R}^n\to\mathbb{R}^n .
$$

For $n=2$ it is often clearer to write the components out,

$$
\frac{\partial x}{\partial t} = f(x, y), \qquad
\frac{\partial y}{\partial t} = g(x, y),
$$

with $\mathbf{u}=(x,y)$. A **fixed point** is again a state where the flow vanishes,
$\mathbf{F}(\mathbf{u}^*) = \mathbf{0}$, i.e. $f(x^*,y^*)=g(x^*,y^*)=0$ simultaneously. As before,
*a solution that starts at a fixed point stays there*, and *a solution that does not start at one
can never reach one*.

The definitions of stable and unstable are word-for-word the same as in one dimension. What is
genuinely new is that **deciding** stability is harder: in 1-D a single number $f'(x^*)$ did the
job; in $n$ dimensions we need a *matrix*.

### 3.1 Two pictures: the vector field and the nullclines

Two visual tools carry almost all the intuition in the plane.

- The **vector field** (or *flow*) attaches the arrow $\mathbf{F}(x,y)=(f,g)$ to each point
  $(x,y)$. Trajectories are the curves that follow these arrows. In Python, `plt.streamplot`
  draws them beautifully.
- The **nullclines** are the curves on which one component of the flow is zero. The
  *$x$-nullcline* is $f(x,y)=0$ (here trajectories move purely vertically); the *$y$-nullcline* is
  $g(x,y)=0$ (trajectories move purely horizontally). **Fixed points are exactly the
  intersections of the two nullclines** — the places where both components vanish at once. Drawing
  nullclines is the fastest way to locate equilibria by eye.

We will use both in every worked example below.

### 3.2 Linear systems and eigenvalues

The cleanest case is a **linear** system,

$$
\frac{d\mathbf{u}}{dt} = A\mathbf{u}, \qquad \mathbf{u}(0)=\mathbf{u}_0,
$$

with $A$ an $n\times n$ matrix. The origin $\mathbf{u}^*=\mathbf{0}$ is always a fixed point, and
(if $A$ is invertible) the only one. Its stability is governed entirely by the **eigenvalues** of
$A$ — the numbers $\lambda$ for which $A\mathbf{v}=\lambda\mathbf{v}$ for some eigenvector
$\mathbf{v}$. Eigenvalues can be real or come in complex-conjugate pairs $\lambda=\alpha\pm i\beta$,
and are found from

$$
\det(A - \lambda I) = 0 .
$$

The single most important theorem of the whole chapter is:

:::{important}
**Stability of a linear system.** The fixed point at the origin of $\dot{\mathbf u}=A\mathbf u$ is
**stable** if *all* eigenvalues of $A$ have **negative real part**, and **unstable** if *any*
eigenvalue has **positive real part**.
:::

The *real part* sets growth or decay; the *imaginary part* sets rotation. This gives the **five
canonical phase portraits** in two dimensions ({numref}`fig-portraits`):

```{figure} figures/linear_portraits.png
:name: fig-portraits

The five behaviours of a 2-D linear system, with eigenvalues in the titles. **Real, both
negative** → *stable node* (everything decays straight in). **Real, both positive** → *unstable
node* (everything grows out). **Real, opposite signs** → *saddle* (attracts along one direction,
repels along another). **Complex with negative real part** → *stable spiral / focus* (decaying
oscillation). **Complex with positive real part** → *unstable spiral* (growing oscillation).
```

Complex eigenvalues always produce **spirals**, and a spiral in the $(u_1,u_2)$ plane corresponds
to **oscillations** in each component as a function of time — the rotation rate is the imaginary
part $\beta$. This is why oscillations and complex eigenvalues are, in practice, the same
phenomenon.

```python
import numpy as np

A = np.array([[-1.0, -10.0],
              [10.0, -1.0]])
eigvals, eigvecs = np.linalg.eig(A)
print(eigvals)          # -1+10j, -1-10j  -> stable spiral
```

### 3.3 Classifying a 2-D fixed point with the trace and determinant

In two dimensions you rarely need to compute the eigenvalues explicitly. Two scalars built
directly from the matrix do everything, because for a $2\times2$ matrix

$$
\det(A) = \lambda_1\lambda_2, \qquad \operatorname{Tr}(A)=\lambda_1+\lambda_2 ,
$$

and the eigenvalues themselves are

$$
\lambda_{1,2} = \frac{T \pm \sqrt{T^2 - 4D}}{2}, \qquad T=\operatorname{Tr}(A),\; D=\det(A).
$$

From these three formulas the entire classification follows by inspection:

| Condition | Type | Stability |
|---|---|---|
| $D < 0$ | **saddle** (eigenvalues real, opposite sign) | unstable |
| $D > 0,\; T < 0,\; T^2 > 4D$ | stable **node** | stable |
| $D > 0,\; T > 0,\; T^2 > 4D$ | unstable **node** | unstable |
| $D > 0,\; T < 0,\; T^2 < 4D$ | stable **spiral** | stable |
| $D > 0,\; T > 0,\; T^2 < 4D$ | unstable **spiral** | unstable |

In words: **the determinant decides node-vs-saddle and the trace decides stable-vs-unstable,
while the discriminant $T^2-4D$ decides node-vs-spiral.** This is summarised in the
**trace–determinant plane** ({numref}`fig-trace-det`), a single map of every possible
two-dimensional linear behaviour, and the picture to keep in your head for the rest of the
chapter.

```{figure} figures/trace_determinant.png
:name: fig-trace-det
:width: 75%

The trace–determinant plane. The horizontal axis is $T=\operatorname{Tr}(A)$, the vertical axis
$D=\det(A)$. The line $D=0$, the line $T=0$, and the parabola $D=T^2/4$ divide the plane into five
regions, one for each portrait. Crossing the parabola turns a node into a spiral; crossing
$T=0$ with $D>0$ is a **Hopf bifurcation**; crossing $D=0$ is a **saddle-node bifurcation**.
```

### 3.4 Nonlinear systems: the Jacobian

Real models are nonlinear, but near a fixed point a smooth nonlinear system *looks linear*. The
matrix that captures this local linearisation is the **Jacobian** — the multidimensional
generalisation of the derivative $f'(x^*)$. For a 2-D system its entries are the partial
derivatives of the flow, evaluated **at the fixed point**:

$$
J \big|_{(x^*,y^*)}=
\begin{bmatrix}
\dfrac{\partial f}{\partial x} & \dfrac{\partial f}{\partial y}\\[2ex]
\dfrac{\partial g}{\partial x} & \dfrac{\partial g}{\partial y}
\end{bmatrix} .
$$

The key theorem says that the local picture is governed by the Jacobian's eigenvalues exactly as
in the linear case:

:::{important}
**Linearisation theorem.** Near a fixed point $\mathbf u^*$, solutions of the nonlinear system
$\dot{\mathbf u}=\mathbf F(\mathbf u)$ behave like solutions of the linear system
$\dot{\mathbf u}=J\mathbf u$ near the origin. In particular: if **all** eigenvalues of $J$ have
negative real part, $\mathbf u^*$ is **stable**; if **any** has positive real part, it is
**unstable**. Complex eigenvalues mean the approach (or escape) is a **spiral**.
:::

So the whole procedure for analysing a nonlinear model is:

1. Find the fixed points by solving $\mathbf F(\mathbf u^*)=\mathbf 0$ (often where the nullclines
   cross).
2. Compute the Jacobian $J$ and evaluate it at each fixed point.
3. Read off the type and stability from $\operatorname{Tr}(J)$, $\det(J)$, and $T^2-4D$.

### 3.5 A reusable analysis toolkit

The following self-contained functions implement steps 1–3 once and for all. Every worked example
below reuses them, changing only the model-specific flow and Jacobian.

```python
import numpy as np
import scipy.integrate
import scipy.optimize

def integrate(flow, y0, t):
    """Numerically integrate dy/dt = flow(y, t) from y0 over time array t."""
    return scipy.integrate.odeint(flow, y0, t)

def find_equilibrium(flow, guess):
    """Solve flow(y) = 0 starting from `guess`; return the root or NaNs."""
    sol, info, ok, msg = scipy.optimize.fsolve(flow, guess, full_output=1)
    return sol if ok == 1 else np.full_like(np.asarray(guess, float), np.nan)

def unique_equilibria(flow, guesses, tol=1e-4):
    """Collect distinct equilibria found from a list of starting guesses."""
    eqs = []
    for g in guesses:
        r = find_equilibrium(flow, g)
        if not np.any(np.isnan(r)) and not any(np.allclose(r, e, atol=tol) for e in eqs):
            eqs.append(r)
    return eqs

def classify(J):
    """Return (name, is_stable) for a 2x2 Jacobian via trace & determinant."""
    T, D = np.trace(J), np.linalg.det(J)
    if D < 0:
        return "saddle", False
    stable = T < 0
    shape = "spiral" if T**2 - 4 * D < 0 else "node"
    return f"{'stable' if stable else 'unstable'} {shape}", stable
```

A natural way to *find* the equilibria of an unfamiliar model is to integrate from many initial
conditions, let the trajectories settle, and use their endpoints as guesses for `fsolve`
(stable equilibria are found this way automatically; unstable ones need guesses placed near them,
e.g. on a nullcline).

:::{admonition} Exercises (systems)
:class: seealso
1. Invent a $2\times2$ matrix, compute its eigenvalues by hand, classify the origin, then check
   with `np.linalg.eig` and by simulating $\dot{\mathbf u}=A\mathbf u$ from a few initial
   conditions.
2. Redo Exercise 1 using only $\operatorname{Tr}(A)$ and $\det(A)$ — no eigenvalues.
3. For the nonlinear system $\dot x = -(1-y^2)x - y,\; \dot y = x$, find the unique fixed point,
   compute the Jacobian there, and classify it. Simulate to confirm.
:::

### 3.6 Bifurcations in systems: saddle-node, Hopf, and limit cycles

Systems inherit the **saddle-node** bifurcation from one dimension: as a parameter changes, two
equilibria can collide and vanish. On the trace–determinant plane this is the moment $\det(J)$
passes through zero ($D=0$), switching a stable node into a saddle.

Systems can also do something *impossible in one dimension*. When a **spiral** changes stability —
when $\operatorname{Tr}(J)$ crosses zero with $\det(J)>0$, so a pair of complex eigenvalues
crosses the imaginary axis — we have a **Hopf bifurcation**. As a stable spiral loses stability it
typically gives birth to a **limit cycle**: an isolated *closed* trajectory that nearby solutions
spiral onto. A limit cycle is a self-sustained oscillation, and it is the mathematical heart of
*rhythmic firing*, *central pattern generators*, and essentially every biological clock. The Hopf
bifurcation is therefore the single most important bifurcation in computational neuroscience, and
we will watch one create spiking in the next section but one.

---

(sec-toggle)=
## 4. Worked example I — Bistability: the genetic toggle switch

Our first full model is the **genetic toggle switch** of Gardner, Cantor and Collins
{cite}`gardner2000`. Two genes each produce a protein that *represses* the other. The
concentrations $u$ and $v$ obey

$$
\frac{du}{dt} = \frac{\alpha}{1+v^{\beta}} - u, \qquad
\frac{dv}{dt} = \frac{\alpha}{1+u^{\beta}} - v .
$$

Here $\alpha$ is the maximum production rate and $\beta$ the **cooperativity** (how sharply one
protein shuts the other off); the $-u$ and $-v$ terms are simple degradation. Although it was
built for synthetic biology, the toggle switch is the prototype for any system that must *store a
discrete memory*: a neuron or circuit that can sit in either a low or a high state and stay there.
The mathematical signature of memory is **bistability** — two coexisting stable fixed points — and
this model lets us see exactly when it appears.

### 4.1 Nullclines, equilibria, and the phase plane

Setting each derivative to zero gives the two nullclines

$$
\dot u = 0 \;\Leftrightarrow\; u = \frac{\alpha}{1+v^{\beta}}, \qquad
\dot v = 0 \;\Leftrightarrow\; v = \frac{\alpha}{1+u^{\beta}} ,
$$

and the equilibria are their intersections. The number of intersections depends on $\beta$
({numref}`fig-toggle-phase`): for low cooperativity the curves cross **once** (a single stable
state — *monostable*), but as $\beta$ grows the nullclines bend enough to cross **three** times,
producing **two stable states separated by a saddle** (*bistable*). The saddle's stable manifold
is the *threshold* that decides which of the two memories the system falls into.

```{figure} figures/toggle_phase.png
:name: fig-toggle-phase

Phase plane of the toggle switch ($\alpha=1$). Teal and orange curves are the $\dot u=0$ and
$\dot v=0$ nullclines; grey streamlines are the flow. **Left ($\beta=2$):** the nullclines cross
once — a single stable node (blue), so the system is *monostable*. **Right ($\beta=10$):** they
cross three times — two stable nodes (blue) flanking a saddle (purple). The system is a genuine
**bistable switch**, and the saddle marks the boundary between the two basins of attraction.
```

```python
from functools import partial
import numpy as np

def toggle(y, t, alpha, beta):
    u, v = y
    return np.array([alpha / (1 + v**beta) - u,
                     alpha / (1 + u**beta) - v])

def toggle_jacobian(u, v, alpha, beta):
    # J = -[[1, a*b*v^(b-1)/(1+v^b)^2], [a*b*u^(b-1)/(1+u^b)^2, 1]]
    return -np.array([[1, alpha*beta*v**(beta-1) / (1 + v**beta)**2],
                      [alpha*beta*u**(beta-1) / (1 + u**beta)**2, 1]])

params = {"alpha": 1.0, "beta": 10.0}
flow = lambda y: toggle(y, 0, **params)

# locate equilibria by relaxing from a spread of initial conditions
guesses = [integrate(partial(toggle, **params), ic, np.linspace(0, 50, 500))[-1]
           for ic in [(.1, 1), (2, 2), (1, 1.3), (2, 3), (2, 1), (1, 2)]]
for eq in unique_equilibria(flow, guesses):
    name, _ = classify(toggle_jacobian(*eq, **params))
    print(f"{name:14s} at (u, v) = ({eq[0]:.3f}, {eq[1]:.3f})")
```

If you prefer to avoid differentiating the Jacobian by hand, `sympy` will do it for you and hand
back a fast numerical function:

```python
import sympy
u, v, alpha, beta = sympy.symbols("u v alpha beta")
F = sympy.Matrix([alpha/(1 + v**beta) - u, alpha/(1 + u**beta) - v])
J = F.jacobian(sympy.Matrix([u, v]))
toggle_jacobian = sympy.lambdify((u, v, alpha, beta), J, dummify=False)
```

### 4.2 The bifurcation diagram and how symmetry shapes it

To see *when* bistability switches on, we track the equilibria as $\beta$ varies. Re-solving from
scratch at every $\beta$ is wasteful; instead we use **natural-parameter continuation** — take the
solution at one $\beta$ as the initial guess for the next. A few lines of code trace each branch:

```python
def continuation(flow_of_param, u0, param_values):
    """Follow a root of flow_of_param(u, p)=0 as p sweeps through param_values."""
    eqs, u = [], np.array(u0, float)
    for p in param_values:
        u = find_equilibrium(lambda y: flow_of_param(y, p), u)
        eqs.append(u.copy())
    return np.array(eqs)

beta_space = np.linspace(0.5, 10, 800)
for start in [(0.5, 0.99), (0.84, 0.84), (0.99, 0.5)]:
    branch = continuation(lambda y, b: toggle(y, 0, 1.0, b), start, beta_space)
    # colour each point by stability via toggle_jacobian + classify ...
```

The result is {numref}`fig-toggle-bif`: a single stable branch for small $\beta$ that **splits**
into two stable branches plus an unstable middle branch as $\beta$ crosses a critical value. Because
the model is symmetric under swapping $u\leftrightarrow v$, this is a **pitchfork bifurcation**.

```{figure} figures/toggle_bifurcation.png
:name: fig-toggle-bif
:width: 78%

Bifurcation diagram of the symmetric toggle switch: the steady-state value of $u$ versus
cooperativity $\beta$. Below $\beta\approx4$ there is one stable state; above it the branch forks
into two stable states (blue) with an unstable saddle (purple) in between. This fork is a
**supercritical pitchfork** — the onset of memory.
```

The pitchfork is a consequence of the model's perfect symmetry, and perfect symmetry is fragile.
If the two genes repress with *different* cooperativities $\beta_1\neq\beta_2$, the pitchfork
**unfolds** into a pair of ordinary saddle-node bifurcations: the fork breaks into one smooth
branch plus a detached fold. Sweeping both $\beta_1$ and $\beta_2$ reveals that the two
saddle-nodes meet at a single point — a **cusp**, the simplest *codimension-2* bifurcation, where
the boundary between "monostable" and "bistable" regions comes to a sharp tip. (The asymmetric
model and its cusp are worked out in full in the accompanying `bistable_systems.ipynb` notebook.)

:::{admonition} Exercises (toggle switch)
:class: seealso
1. Reproduce {numref}`fig-toggle-phase` and confirm by counting nullcline intersections that
   $\beta=2$ is monostable while $\beta=10$ is bistable. Roughly where does the transition happen?
2. Starting from two nearby initial conditions on opposite sides of the saddle, simulate and show
   that they converge to *different* stable states. This is the switch *remembering* its input.
3. Break the symmetry: set $\beta_1=6$, sweep $\beta_2$, and identify the saddle-node bifurcations
   in the resulting diagram.
:::

---

(sec-fhn)=
## 5. Worked example II — Excitability: the FitzHugh–Nagumo neuron

The **FitzHugh–Nagumo (FHN)** model is a two-variable caricature of the Hodgkin–Huxley equations
that nonetheless reproduces the essential repertoire of a neuron — rest, threshold, spike, and
repetitive firing. It reads

$$
\frac{dv}{dt} = v - v^3 - w + I_\text{ext}, \qquad
\tau\,\frac{dw}{dt} = v - a - b\,w ,
$$

where $v$ is a (dimensionless) membrane potential, $w$ a slow **recovery** variable, $\tau$ its
(large) time constant, and $I_\text{ext}$ an injected current. The cubic in $v$ supplies the fast
self-amplifying upstroke of a spike; the slow $w$, pulled along by $v$ and decaying as $-bw$,
provides the recovery that drags the cell back toward rest. Three experimental facts motivate the
structure: the cell rests quietly; a *small* perturbation decays back to rest; but a perturbation
past a **threshold** triggers a large excursion — a spike — before returning. A system with this
"ignore small kicks, amplify big ones" property is called **excitable**.

### 5.1 Nullclines, equilibria, and the phase plane

The nullclines are

$$
\dot v = 0 \;\Leftrightarrow\; w = v - v^3 + I_\text{ext}\quad\text{(cubic)}, \qquad
\dot w = 0 \;\Leftrightarrow\; w = \tfrac{1}{b}(v - a)\quad\text{(line)} .
$$

Equilibria satisfy a single cubic in $v$, which we can solve exactly with `numpy.roots`:

$$
v^3 + v\!\left(\tfrac{1}{b}-1\right) - \left(\tfrac{a}{b}+I_\text{ext}\right) = 0,
\qquad w^* = v^* - (v^*)^3 + I_\text{ext} .
$$

```python
import numpy as np
from functools import partial

def fitzhugh_nagumo(x, t, a, b, tau, I):
    v, w = x
    return np.array([v - v**3 - w + I,
                     (v - a - b * w) / tau])

def fhn_equilibria(a, b, tau, I):
    coeffs = [1, 0, 1/b - 1, -a/b - I]               # v^3 + (1/b - 1) v - (a/b + I)
    return [[r.real, r.real - r.real**3 + I]
            for r in np.roots(coeffs) if abs(r.imag) < 1e-9]

def fhn_jacobian(v, w, a, b, tau, I):
    return np.array([[1 - 3*v**2, -1.0],
                     [1/tau,      -b/tau]])

p = {"a": -0.3, "b": 1.0, "tau": 20, "I": 0.0}
for v, w in fhn_equilibria(**p):
    name, _ = classify(fhn_jacobian(v, w, **p))
    print(f"{name:14s} at v* = {v:+.3f}")
```

{numref}`fig-fhn-phase` shows the phase plane in two regimes. **At rest ($I=0$)** the line crosses
the cubic once, at a single **stable** fixed point: small displacements relax back along the
streamlines, but a displacement past the middle branch of the cubic launches a long loop — a
single spike — before returning. **With enough current ($I=0.4$)** the intersection moves onto the
unstable middle branch of the cubic, the fixed point becomes **unstable**, and the trajectory can
no longer settle: it is captured by a **limit cycle** (red), and the neuron fires *repetitively*.

```{figure} figures/fhn_phase.png
:name: fig-fhn-phase

FitzHugh–Nagumo phase plane ($a=-0.3,\,b=1.0,\,\tau=20$). The teal cubic is the $\dot v=0$
nullcline, the orange line the $\dot w=0$ nullcline; their intersection is the fixed point.
**Left ($I=0$):** a single *stable* fixed point (blue) — the resting state; trajectories from
small kicks (grey) decay back. **Right ($I=0.4$):** the fixed point is now *unstable* (red) and is
encircled by a *stable limit cycle* (red loop) — the neuron spikes periodically.
```

The same two regimes seen as voltage traces are unmistakable
({numref}`fig-fhn-time`): a sub-threshold bump that decays on the left, a regular spike train on
the right.

```{figure} figures/fhn_timeseries.png
:name: fig-fhn-time

Membrane potential $v(t)$ for the same two cases. **Left ($I=0$):** small kicks decay; a larger
kick produces one spike and then quiescence (excitable but not oscillating). **Right ($I=0.4$):**
sustained, periodic spiking — the time-domain face of the limit cycle.
```

### 5.2 What the current does: a Hopf bifurcation and excitation block

Sweeping $I_\text{ext}$ and classifying the fixed point at each value produces the bifurcation
diagram in {numref}`fig-fhn-bif`. As $I$ increases, the resting fixed point loses stability at a
**Hopf bifurcation** ($\operatorname{Tr}(J)$ crosses zero), the limit cycle is born, and the neuron
begins to fire; at a *second* Hopf the fixed point regains stability and the firing stops. This
upper transition is **excitation block** (or *depolarisation block*): inject too much current and a
real neuron, counter-intuitively, goes *silent*, pinned at a high, depolarised steady state.
Spiking lives only in the band between the two Hopf points.

```{figure} figures/fhn_bifurcation.png
:name: fig-fhn-bif
:width: 78%

FitzHugh–Nagumo bifurcation diagram: steady-state potential $v^*$ versus injected current
$I_\text{ext}$ ($b=1.0$). Blue marks a stable fixed point, red an unstable one. The fixed point is
stable at low and high current and **unstable in the shaded band**, bounded by two Hopf
bifurcations. Inside this band the stable rest state is replaced by a limit cycle — this is the
neuron's firing range. Past the upper Hopf the cell falls silent again (excitation block).
```

```python
import matplotlib.pyplot as plt

I_values = np.linspace(-0.2, 1.0, 500)
for I in I_values:
    p = {"a": -0.3, "b": 1.0, "tau": 20, "I": I}
    for v, w in fhn_equilibria(**p):
        name, stable = classify(fhn_jacobian(v, w, **p))
        plt.plot(I, v, ".", color="C0" if stable else "C3", ms=4)
plt.xlabel(r"$I_{\rm ext}$"); plt.ylabel(r"$v^*$"); plt.show()
```

A two-parameter sweep over $(I, b)$ refines this picture, carving the plane into *monostable*,
*bistable*, and *periodic-firing* regions — a compact map of every behaviour the model can
produce. That analysis, again, is carried out in the companion `excitable_systems.ipynb` notebook.

(sec-nonauto)=
### 5.3 Beyond autonomy: time-varying input and noise

Holding $I_\text{ext}$ constant kept the system autonomous and let us use every tool above. Two
realistic extensions break that assumption but build directly on the same model.

- **Non-autonomous forcing.** Let the current depend on time, $I_\text{ext}=I_\text{ext}(t)$ — a
  step, a sinusoid, a noisy waveform. The system is no longer autonomous, so "fixed points" must
  be reinterpreted, but the phase-plane intuition still guides what the trajectory does as the
  nullclines slide up and down with the input.
- **Stochastic dynamics.** Replacing the ODE with a **stochastic differential equation (SDE)**,
  $dY_t = f(Y_t,t)\,dt + g(Y_t,t)\,dB_t$, adds Brownian noise $B_t$. Numerically this is handled by
  the **Euler–Maruyama** scheme — forward Euler with an extra noise increment scaled by
  $\sqrt{dt}$:

```python
def euler_maruyama(drift, diffusion, y0, t):
    """Integrate dY = drift(Y,t) dt + diffusion(Y,t) dB, B Brownian."""
    y = np.zeros((len(t), len(y0)))
    y[0] = y0
    for n, dt in enumerate(np.diff(t), 1):
        noise = np.random.normal(0.0, np.sqrt(dt))
        y[n] = y[n-1] + drift(y[n-1], t[n]) * dt + diffusion(y[n-1], t[n]) * noise
    return y
```

  Near a stable rest state just below the firing threshold, noise produces the hallmark of
  excitability: most fluctuations decay, but the occasional large excursion crosses threshold and
  fires a full spike — *noise-induced spiking*, a phenomenon you can explore with the code above.

:::{admonition} Exercises (FitzHugh–Nagumo)
:class: seealso
1. Reproduce {numref}`fig-fhn-phase`. Starting just below and just above the middle branch of the
   cubic nullcline, show that one initial condition decays while the other produces a full spike —
   demonstrating the *threshold*.
2. Locate the lower Hopf bifurcation in $I$ (for $b=1.0$) by finding where $\operatorname{Tr}(J)$
   changes sign at the fixed point. Confirm that spiking begins there.
3. Increase $I$ past the upper Hopf and show that the cell stops firing (excitation block).
4. Drive the model with the Euler–Maruyama code and a small noise amplitude near rest. Estimate
   how the spike rate depends on the noise level.
:::

---

## 6. Summary

| Question | One dimension | Systems ($n\ge2$) |
|---|---|---|
| What is a fixed point? | $f(x^*)=0$ | $\mathbf F(\mathbf u^*)=\mathbf 0$ (nullclines cross) |
| How do I picture the dynamics? | phase line | vector field + nullclines (phase plane) |
| What decides stability? | sign of $f'(x^*)$ | real parts of the eigenvalues of the **Jacobian** $J$ |
| Quick 2-D test | — | $\det J$: node/saddle; $\operatorname{Tr}J$: stable/unstable; $T^2-4D$: node/spiral |
| Typical bifurcations | saddle-node | saddle-node **and Hopf** (creates limit cycles) |
| Neuroscience meaning | threshold / switching on | bistable memory; rhythmic / repetitive firing |

The single most transferable idea is that **a model's behaviour is organised by its fixed points
and the bifurcations that create, destroy, or destabilise them.** Locate the equilibria, linearise
to get the Jacobian, read off stability from the eigenvalues, and then ask how that structure
reorganises as a biological parameter — a current, a conductance, a coupling strength — is turned
up or down. Bistability (two stable states) underlies discrete memory; a Hopf bifurcation and its
limit cycle underlie every kind of neural oscillation. The rest of this book returns to these two
motifs again and again.

---

## References

```{bibliography}
:filter: docname in docnames
```

<!-- If your book does not yet use sphinxcontrib-bibtex, replace the directive above with a
plain list and add these entries to your references.bib:

@book{strogatz, author={Strogatz, Steven H.}, title={Nonlinear Dynamics and Chaos},
  publisher={Westview Press}, year={2015}, edition={2nd}}
@book{rosenbaum2022, author={Rosenbaum, Robert},
  title={Modeling Neural Circuits Made Simple with Python}, publisher={MIT Press}, year={2022}}
@article{gardner2000, author={Gardner, T. S. and Cantor, C. R. and Collins, J. J.},
  title={Construction of a genetic toggle switch in Escherichia coli},
  journal={Nature}, volume={403}, pages={339--342}, year={2000}}
-->
