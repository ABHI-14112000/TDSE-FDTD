# TDSE-FDTD (MATLAB)

This repository contains a MATLAB example for solving the **time-dependent Schrödinger equation (TDSE)** for a **1D quantum dot** using a finite-difference grid (FDTD-style discretization in space) with **Crank-Nicolson** time integration.

## File

- `tdse_quantum_dot_fdtd.m`: complete simulation script.

## What it models

- Electron in an effective-mass approximation.
- 1D finite potential well (quantum dot) in a finite simulation box.
- Initial Gaussian wave packet inside the dot.
- Time evolution with stable implicit update.

## Run

Open MATLAB and run:

```matlab
tdse_quantum_dot_fdtd
```

The script will:

1. Build the Hamiltonian matrix from second-order finite differences.
2. Evolve the wavefunction in time.
3. Plot probability density and wavefunction components.
4. Print final normalization and total energy.

## Notes

- You can change `wellWidth`, `V0_eV`, `dt`, `Nx`, and `Nt` to study accuracy and physical behavior.
- For stronger confinement or sharper features, increase `Nx` and reduce `dt`.
