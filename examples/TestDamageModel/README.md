# Smoothed Continuum Damage Model — Single-Element Test

This example exercises the `SmoothDamageModel` (`fingal.damage`, the smoothed
continuum damage model of Mondal, Olsen-Kettle & Gross, *Computers and
Geotechnics* **122** (2020) 103505) on a single finite element before it is used
on real meshes.

The driver [`test1.py`](./test1.py) performs two checks:

1. **Constitutive check** — the damage law `D(kappa)` (Eq. 8) is evaluated over
   `[kappa0, kappa_c]` and plotted (cf. Fig. 2 of the paper).
2. **FEM check** — a single tri-linear hexahedral element is loaded in
   displacement-controlled uniaxial compression with free lateral faces (minimal
   corner pins allow Poisson expansion). A staggered scheme updates the damage
   and re-solves equilibrium, producing the axial stress–strain and
   damage–strain response.

On a single element the non-local Helmholtz smoothing is the identity (uniform
field, Neumann boundary), so this exercises the full local model and the
damage/elasticity coupling.

### Run

    PYTHONPATH=../../bin run-escript test1.py

Outputs:

- `damage_vs_strain.png` — the constitutive damage curve.
- `single_element_response.png` — stress–strain and damage evolution.

The run also prints the per-step response and, at the end, the peak stress and
final damage (≈ 30.4 MPa and `D = 1.0` for the shipped parameters).

### Loading and adaptive step-size control

The loading path is driven by `model.runLoading(set_bc, nsteps, ...)`. It is
parameterised by a **load factor** in `[0, 1]` with nominal increment
`1/nsteps`; the load level at factor `f` is prescribed through
`set_bc(pde, f*nsteps, nsteps)`, so at the nominal factors `k/nsteps` this is
just the plain `set_bc(pde, k, nsteps)` that the example defines. Because
`set_bc` is also evaluated at intermediate factors, it must accept a *fractional*
`step` argument — the standard `u_bc * step / nsteps` form (used in
[`test1.py`](./test1.py)) does.

Each (sub-)step calls `solveLoadStep`, which runs the staggered
equilibrium/damage iteration to the tolerance `tol` (default `1e-6`) with at most
`max_iter` iterations (default `20`).

`runLoading` accepts the following step-control options:

| Option | Default | Meaning |
| --- | --- | --- |
| `nsteps` | — | number of nominal load increments (nominal increment `1/nsteps`). |
| `tol` | `1e-6` | convergence tolerance on the damage change per staggered iteration. |
| `max_iter` | `20` | maximum staggered iterations per (sub-)step. |
| `adaptive` | `True` | enable automatic bisection sub-stepping (see below). |
| `min_increment` | `1e-3/nsteps` | smallest allowed load-factor increment; reaching it without convergence terminates the run. |

**Adaptive bisection (`adaptive=True`, the default).** When a (sub-)step does not
converge within `max_iter` staggered iterations, the committed state
(`self.u`, `self.stress`, `self.kappa`, `self.D`) is rolled back, the load
increment is **halved** and the (sub-)step is retried. After a successful step
the increment is grown back toward the nominal `1/nsteps`. If a (sub-)step still
does not converge once the increment has been reduced to `min_increment`, the
calculation is **terminated** with a `RuntimeError` (rather than silently
accepting an unconverged step). Bisection relies on the incremental-residual
formulation of `solveLoadStep`, which carries the stress state between steps.

Rollback requires that only the load level changes between (sub-)steps, so
`set_bc` must be a pure function of the load factor (no external state).

**Fixed schedule (`adaptive=False`).** The classic fixed `nsteps` schedule is
used with no sub-stepping; a non-converged step is committed anyway and
`solveLoadStep` logs a warning. This reproduces the original behaviour.

Examples:

    # default: adaptive on, terminates if it cannot converge at 1e-3/nsteps
    model.runLoading(set_bc, nsteps=80, callback=record)

    # allow finer sub-stepping before giving up
    model.runLoading(set_bc, nsteps=80, callback=record, min_increment=1e-6)

    # disable adaptivity (fixed 80-step schedule, warn on non-convergence)
    model.runLoading(set_bc, nsteps=80, callback=record, adaptive=False)

Progress and step-control messages (bisection, non-convergence, termination) are
emitted through the `fingal.SmoothDamageModel` logger; enable them with

    import logging
    logging.basicConfig(level=logging.INFO)

The convergence flag and final damage change of the most recent step are also
available as `model.converged` and `model.last_change`.
