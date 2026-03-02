# TDSE-FDTD (MATLAB)

MATLAB example for solving the **1D time-dependent Schrödinger equation (TDSE)** for a **quantum dot** using an **explicit FDTD leapfrog method**.

## File

- `tdse_quantum_dot_fdtd.m` — self-contained script.

## Model

- Effective-mass electron in 1D.
- Finite square-well quantum dot inside a barrier.
- Real/imaginary split of wavefunction with staggered-time FDTD update.
- Optional absorbing boundary layer (CAP-like damping) to reduce edge reflections.

## Run

In MATLAB:

```matlab
tdse_quantum_dot_fdtd
```

## Key parameters to tune

- Spatial grid: `Nx`, `L`
- Time stepping: `safety`, `Nt`
- Dot potential: `wellWidth`, `V0eV`
- Initial packet: `x0`, `sigma`, `k0`
- Absorber: `absWidth`, `etaMax`

## Notes

- The script computes and prints a stability reference time step (`dtStability`) and uses `dt = safety*dtStability`.
- If you see numerical noise, reduce `safety` or increase `Nx` with a correspondingly smaller `dt`.
