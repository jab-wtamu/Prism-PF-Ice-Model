# Demange Ice Growth Model — Fern Dendrite (PRISMS-PF 4.0.0)

This application provides a 3D Fern dendrite simulation on a
**320 × 320 × 51.2** physical domain with **adaptive mesh refinement (AMR) enabled**.
The configuration in `configs/paper_fern_dendrite_amr_t500.prm` specifies the
model parameters and numerical setup.

## Files

- `main.cc` — declares the model fields and computational workflow.
- `custom_pde.h` — corrected Demange phase-field model and initial conditions.
- `demange_anisotropy.h` — anisotropy calculations.
- `CMakeLists.txt` — build definition.
- `configs/paper_fern_dendrite_amr_t500.prm` — the **single** Fern configuration.
- `build.sh` — portable *convenience* build wrapper for an installed PRISMS-PF/deal.II environment.
- `run_fern.sb` — portable *convenience* ACES-style Slurm launcher.
- `.gitignore` — prevents accidental commits of build and run outputs.

The build and Slurm scripts are example wrappers. Configure the environment,
library paths, and scheduler options for the target system before use.

## Numerical setup

| Setting | Value |
|---|---|
| Geometry | 320 × 320 × 51.2 |
| Initial mesh subdivisions | 25 × 25 × 4 |
| Initial global refinement | 3 |
| Spatial method | Q2 finite elements in 3D |
| AMR | ON, levels 0–3; remeshing every 2,500 steps |
| Refinement indicator | Phase-field `phi` in [-0.99, 0.99] |
| Time step | 0.0004 |
| Configured final time | T=500 (1,250,000 time steps) |
| Boundary conditions | Natural (zero normal flux) on all outer faces |
| Initial seed | Disc radius 6.4; half-thickness 1.6 |

**Relation to the published model:** the 320 × 320 × 51.2 *physical domain* corresponds to the
published Fern setup. However, this implementation uses adaptive Q2 finite
elements, explicit time stepping, and natural boundaries rather than the
paper's fixed Fourier-spectral discretization, semi-implicit stepping, and
periodic boundaries. The choices for `half_convention`, `q_normalization`, and anisotropy regularization
are explicit model settings; this is **not an exact numerical reproduction**
of the published simulation.

The model target of T=500 is the configured stopping time; it does **not**
mean the original AMR-ON production run completed T=500 within its 72-hour
walltime limit.

## Build and launch

Build with a compatible PRISMS-PF 4.0.0 and deal.II installation:

```bash
export PRISMS_PF_INSTALL=/path/to/PRISMS-PF/install
export DEAL_II_INSTALL=/path/to/deal.II/install
bash build.sh
```

From the application folder, submit a simulation with the appropriate Slurm account:

```bash
sbatch --account=YOUR_ACES_ACCOUNT run_fern.sb configs/paper_fern_dendrite_amr_t500.prm Fern_320_AMR_ON
```

Outputs from the example launcher go under `runs/Fern_320_AMR_ON/`.
The launcher does not overwrite a preexisting run directory.
