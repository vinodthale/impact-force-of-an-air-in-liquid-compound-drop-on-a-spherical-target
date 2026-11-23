# Impact force of an air-in-liquid compound drop on a spherical target

Basilisk C code for an axisymmetric water drop containing an air bubble and impacting a rigid spherical target. The simulation saves the flow field. A separate extractor integrates pressure and viscous tractions over the target to recover the axial impact force.

## Supplied case

| Quantity | Value |
| --- | --- |
| Drop diameter, $D$ | 2.050 mm |
| Impact speed, $U_0$ | 0.50 m/s |
| Bubble diameter, $D_b/D$ | 0.50 |
| Initial geometry | Spherical drop and concentric spherical bubble |
| Target | Hemisphere with radius $0.629960524947428D$ |
| Wall model | No slip with the ideal non-wetting gas-side stencil |
| Initial surface gap | $0.04D$ |
| Meridional domain width | $4D$ |
| Initial and maximum refinement levels | 9 and 11 |
| Snapshot interval | $0.01D/U_0$ |

The case inputs are in `constants.h`. The 90° angle describes the spherical-cap geometry. The finite contact-angle setting is inactive in the supplied gas-cushion configuration.

## Run the simulation

Use Linux, Bash, Basilisk C with `qcc`, and an MPI C compiler and launcher. The force-extraction script reads the velocity and bubble-size labels from the case directory name. Clone the supplied case into the following layout and run it with the local MPI launcher (six ranks by default).

```bash
mkdir -p V0.50
git clone https://github.com/vinodthale/impact-force-of-an-air-in-liquid-compound-drop-on-a-spherical-target.git V0.50/DB0.50
cd V0.50/DB0.50
bash LaptopMPIV0.sh
```

On a Slurm cluster, compile with `bash compile_Basilisk.sh` and submit `sbatch ClusterMPI.sh`. The supplied job requests 64 MPI ranks. Adjust the partition, MPI launch option, and module names to match your site. The scripts can use an already configured Basilisk and MPI environment. The alternative cluster setup reads `$BASILISK_ROOT/scripts/env.sh` and expects the Basilisk source in `$BASILISK_ROOT/app/src`.

## Extract the axial impact force

From the same case directory, run the following command after the simulation finishes. The completeness check verifies the saved time range and checks for gaps between snapshots.

```bash
bash RJetandpressure_embded.sh --require-complete
```

The extractor reads `intermediate/snapshot-*` and writes `Post/impact_force_Vel0.50_DB0.50_90_deg.txt`. The four whitespace-delimited columns are the raw nondimensional simulation time, axial pressure force, axial viscous force, and their sum. Time is in units of $D/U_0$ and force is in units of $\rho_l U_0^2 D^2$. Positive force is in the impact direction. The time origin is the initial drop position. Nominal contact occurs at raw time 0.04 for the supplied case.

## Source files

| File | Purpose |
| --- | --- |
| `Bdropimpact.c`, `constants.h` | Impact simulation and physical and numerical inputs |
| `Jetandpressure_embded.c`, `axi_embed_impact_force.h` | Snapshot restoration and integration of axial traction over the embedded target |
| `BpostBview.c`, `Bview.sh` | Optional Basilisk View rendering of saved interfaces |
| `stage_b_*`, `huang_axi_locked_include/`, `huang_original/` | Geometry and contact-line routines for the finite-angle configuration |
