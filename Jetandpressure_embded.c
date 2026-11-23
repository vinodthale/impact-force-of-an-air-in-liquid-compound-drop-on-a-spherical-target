/* Offline axial force on the curved target, from existing Basilisk dumps.
   The case launcher supplies the saved filenames in chronological order. */
#define WRITE_PARAMETERS_FILE 0
#include "constants.h"
#include "axi_embed_impact_force.h"

struct CFDValues cfdbv;

/* Same wall/far-field boundary conditions as the production solver. */
u.n[embed] = dirichlet (0.);
u.t[embed] = dirichlet (0.);
p[embed]   = neumann (0.);

u.n[right] = neumann (0.);
p[right]   = dirichlet (0.);

u.n[top] = neumann (0.);
p[top]   = dirichlet (0.);

pressure[embed] = neumann (0.);
pressure[right] = dirichlet (0.);
pressure[top]   = dirichlet (0.);

/* Rebuild the same fixed spherical-cap EMBED geometry as Bdropimpact.c. */
static void reconstruct_embed_geometry (void)
{
  const double solid_center = cfdbv.rs - cfdbv.xm;

  solid (cs, fs, sq(x - solid_center) + sq(y) - sq(cfdbv.rs));

  fractions_cleanup (cs, fs);

  restriction ({cs, fs});
  boundary ({cs, fs});

  cm_update (cm, cs, fs);
  fm_update (fm, cs, fs);

  restriction ({cm, fm});
  boundary ({cm, fm});
}

/* Restore the saved fields and reconstruct the face geometry before using
   the case's filtered phase fraction to rebuild the face viscosity. */
static bool write_snapshot_force (const char * path, FILE * history)
{
  if (!restore (file = path)) {
    fprintf (stderr, "Cannot restore %s\n", path);
    return false;
  }
  char expected[128];
  snprintf (expected, sizeof(expected), "intermediate/snapshot-%5.4f", t);
  if (!isfinite(t) || strcmp(path, expected) != 0) {
    fprintf (stderr, "%s has inconsistent stored time %.17g\n", path, t);
    return false;
  }

  reconstruct_embed_geometry();
  event ("properties");
  boundary ({pressure, f, sf, u});
  const AxiTopPressureDiagnostic top = axi_top_adjacent_pressure (pressure);
  if (top.invalid_faces || !isfinite(top.mean)) {
    fprintf (stderr, "Invalid far-field pressure in %s\n", path);
    return false;
  }
  const AxiEmbedImpactForce force = axi_embed_impact_force (pressure, u, mu, top.mean);
  if (force.invalid_fragments || force.fragments <= 0 ||
      !isfinite(force.pressure_x) || !isfinite(force.viscous_x)) {
    fprintf (stderr, "Invalid wall traction in %s (%d fragments, %d invalid)\n",
             path, force.fragments, force.invalid_fragments);
    return false;
  }
  const double fp = -force.pressure_x, fmu = -force.viscous_x;
  return fprintf (history, "%.4f %.10f %.10f %.10f\r\n", t, fp, fmu, fp + fmu) >= 0;
}

/* Check the small filename list before restoring potentially large MPI dumps. */
static bool check_snapshot_coverage (FILE * list, double stop, bool require_complete)
{
  char path[128];
  double first = NAN, last = NAN;
  bool gap = false;
  while (fgets(path, sizeof(path), list)) {
    double saved;
    char trailing;
    if (sscanf(path, "intermediate/snapshot-%lf %c", &saved, &trailing) != 1 ||
        !isfinite(saved) || (isfinite(last) && saved <= last)) {
      fprintf (stderr, "Invalid or unordered snapshot list: %s", path);
      return false;
    }
    if (!isfinite(first))
      first = saved;
    else if (saved - last > SAVE_FILE_EVERY + 0.0001) {
      fprintf (stderr, "Gap between saved snapshots at %.4f and %.4f\n", last, saved);
      gap = true;
    }
    last = saved;
  }
  if (ferror(list) || !isfinite(first)) {
    fprintf (stderr, "Cannot read the saved snapshot list.\n");
    return false;
  }
  if (require_complete && (gap || fabs(first) > 0.00005 ||
                           fabs(last - stop) > 0.00005)) {
    fprintf (stderr, "Incomplete snapshot history: expected 0.0000 through %.4f, got %.4f through %.4f.\n",
             stop, first, last);
    return false;
  }
  if (fseek(list, 0L, SEEK_SET)) {
    perror ("Cannot rewind snapshot list");
    return false;
  }
  return true;
}

int main (int argc, char ** argv)
{
  /* Folder labels preserve V0.18 even though its physical input speed
     is 0.18889446348220329 m/s. The force calculation uses constants.h. */
  if ((argc != 3 && argc != 4) ||
      (argc == 4 && strcmp(argv[3], "--require-complete") != 0)) {
    fprintf (stderr, "Usage: extractor VELOCITY_LABEL DB_LABEL [--require-complete]\n");
    return EXIT_FAILURE;
  }
  const bool require_complete = argc == 4;
  char name[128];
  int n = snprintf (name, sizeof(name), "impact_force_Vel%s_DB%s_%d_deg.txt",
                    argv[1], argv[2], (int)(SOLID_SHAPE_ANGLE_DEGREES + 0.5));
  if (n < 0 || (size_t)n >= sizeof(name))
    return EXIT_FAILURE;
  FILE * snapshots = fopen ("snapshots.list", "r");
  if (!snapshots) {
    fprintf (stderr, "Run RJetandpressure_embded.sh to prepare the snapshot list.\n");
    return EXIT_FAILURE;
  }
  /* The launcher runs here in a fresh private directory. Reserve a new
     report before the shared setup opens its report for writing. */
  FILE * new_report = fopen (FILENAME_PARAMETERS, "wx");
  if (!new_report) {
    perror ("Cannot reserve a new, private parameter report");
    fclose (snapshots);
    return EXIT_FAILURE;
  }
  fclose (new_report);
  /* The two folder labels are not solver parameter overrides. */
  numericalmainvalues (argv, 1, &cfdbv);
  if (!check_snapshot_coverage (snapshots, cfdbv.timecontact + cfdbv.timeend,
                                require_complete)) {
    fclose (snapshots);
    return EXIT_FAILURE;
  }
  rho1 = cfdbv.rhoL; rho2 = cfdbv.rhoG;
  mu1  = cfdbv.muL;  mu2  = cfdbv.muG;
  size (cfdbv.domainsize);
  origin (0., 0.);
  init_grid (1 << LEVELmin);
  event ("metric");
  event ("defaults");

  FILE * history = fopen (name, "wx");
  if (!history) {
    perror ("Cannot open force history");
    fclose (snapshots);
    free_grid();
    return EXIT_FAILURE;
  }
  fprintf (history, "# Solid Shape Angle [deg]: %.10f\n", (double) SOLID_SHAPE_ANGLE_DEGREES);
  fputs ("# t_star Fp_star Fmu_star Ftotal_star\n"
         "# Raw time in D/U0 units. Axial pressure, viscous and total impact forces"
         " in rho_l*U0^2*D^2 units. Ftotal_star = Fp_star + Fmu_star.\n", history);

  char path[128];
  int rows = 0;
  double first = 0., last = 0.;
  bool ok = true;
  while (fgets(path, sizeof(path), snapshots)) {
    path[strcspn(path, "\r\n")] = '\0';
    if (!write_snapshot_force (path, history)) {
      ok = false;
      break;
    }
    if (!rows)
      first = t;
    last = t;
    rows++;
  }

  if (ferror(snapshots) || ferror(history) || !rows) {
    fprintf (stderr, "Could not read or write a force history with data rows.\n");
    ok = false;
  }
  fclose (snapshots);
  if (fclose(history) != 0) {
    perror ("Cannot finish force history");
    ok = false;
  }
  free_grid();
  if (!ok)
    return EXIT_FAILURE;
  fprintf (stderr, "Extracted %d snapshots, t*=%.4f through %.4f (configured stop %.4f).\n",
           rows, first, last, cfdbv.timecontact + cfdbv.timeend);
  return EXIT_SUCCESS;
}
