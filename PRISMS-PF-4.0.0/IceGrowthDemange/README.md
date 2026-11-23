# Demange Ice Growth Model — PRISMS-PF 4.0.0

This application models three-dimensional snow-crystal growth using a phase-field
formulation. Two adaptive-mesh-refinement (AMR) configurations are provided:
a Fern dendrite and a Hollow Prism.

## Files

- `main.cc` — declares the model fields and computational workflow.
- `custom_pde.h` — Demange phase-field equations and initial conditions.
- `demange_anisotropy.h` — anisotropy calculations.
- `CMakeLists.txt` — build definition.
- `configs/paper_fern_dendrite_amr_t500.prm` — Fern simulation parameters.
- `configs/hollow_Gamma3_U128_T100_AMR_ON_dt2e4.prm` — Hollow Prism simulation parameters.
- `build.sh` — example build script for PRISMS-PF and deal.II.
- `run_fern.sb` — example Slurm launcher for an ACES-style environment.
- `.gitignore` — prevents accidental commits of build and run outputs.

The build and Slurm scripts are example wrappers. Configure the environment,
library paths, and scheduler options for the target system before use.

## Fern dendrite numerical setup

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

T=500 is the configured stopping time. A scheduler walltime limit may prevent
an individual job from reaching that target.

## Hollow Prism numerical setup

| Setting | Value |
|---|---|
| Geometry | 51.2 × 51.2 × 102.4 |
| Initial mesh subdivisions | 4 × 4 × 8 |
| Initial global refinement | 5 |
| Spatial method | Q2 finite elements in 3D |
| AMR | ON, levels 0–5; remeshing every 2,500 increments |
| Refinement indicator | Phase-field `phi` in [-0.99, 0.99] |
| Time step | 0.0002 |
| Configured final time | T=100 (500,000 increments) |
| Boundary conditions | Natural (zero normal flux) on all outer faces |
| Initial seed | Disc radius 3.2; half-thickness 1.6 |

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

To run the Hollow Prism configuration, use the same application executable
with the corresponding parameter file. The included Slurm launcher retains
the job name `Fern`; adjust its job label when using the Hollow configuration.

Outputs from the Fern example launcher go under `runs/Fern_320_AMR_ON/`.
The launcher does not overwrite a preexisting run directory.
