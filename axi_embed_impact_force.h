/*
 * Axisymmetric traction diagnostic for a two-dimensional Basilisk EMBED
 * surface.  This header is diagnostic-only: it reads p, u, mu, cs, fs and fm
 * and never modifies solver fields.
 *
 * The local Basilisk embed_force() integrates meridional fragment length but
 * does not call metric_embed_factor().  For AXI, each fragment is therefore
 * revolved explicitly using dA = 2*pi*r_b*dl, where r_b is the fragment-
 * barycentre radius and dl is its meridional length.
 */

#ifndef AXI_EMBED_IMPACT_FORCE_H
#define AXI_EMBED_IMPACT_FORCE_H

#if dimension != 2
# error "axi_embed_impact_force.h requires a two-dimensional grid"
#endif
#if !AXI
# error "axi_embed_impact_force.h requires axi.h"
#endif
#if !EMBED
# error "axi_embed_impact_force.h requires embed.h"
#endif

typedef struct {
  double pressure_x;             // Signed fluid-on-solid pressure force.
  double viscous_x;              // Signed fluid-on-solid viscous force.
  double area;                   // Revolved embedded area, integral dA.
  double max_surface_pressure;   // Maximum interpolated p - p_reference.
  int fragments;                 // Number of cut-cell surface fragments.
  int invalid_fragments;         // Non-finite or degenerate fragments.
} AxiEmbedImpactForce;

typedef struct {
  double mean;
  double area;
  int invalid_faces;
} AxiTopPressureDiagnostic;

/**
 * Integrate the axial Cauchy traction over all embedded fragments.
 *
 * This follows the pressure interpolation, face-metric viscosity recovery,
 * embedded Dirichlet gradient and traction signs of Basilisk's embed_force().
 * The only formulation addition is the missing AXI surface measure 2*pi*r_b.
 * The result has Basilisk nondimensional force units rho_l*U0^2*D^2.
 */
static inline AxiEmbedImpactForce
axi_embed_impact_force (scalar pressure_field, vector velocity,
                        face vector viscosity, double p_reference)
{
  double pressure_x = 0., viscous_x = 0., area_sum = 0.;
  double max_surface_pressure = -HUGE;
  /* Double counters avoid a known false-positive -Warray-bounds warning in
   * this checkout's generated MPI reducer for integer foreach reductions. */
  double fragments = 0., invalid_fragments = 0.;

  foreach (reduction(+:pressure_x) reduction(+:viscous_x)
           reduction(+:area_sum) reduction(max:max_surface_pressure)
           reduction(+:fragments) reduction(+:invalid_fragments), nowarning) {
    if (cs[] > 0. && cs[] < 1.) {
      fragments++;

      coord barycentre, normal;
      const double meridional_fraction =
        embed_geometry (point, &barycentre, &normal);
      const double radius = y + barycentre.y*Delta;
      const double dA = 2.*pi*radius*meridional_fraction*Delta;

      if (!(meridional_fraction > 0.) || !(radius >= 0.) || !(dA > 0.) ||
          !isfinite(meridional_fraction) || !isfinite(radius) ||
          !isfinite(dA) || !isfinite(normal.x) || !isfinite(normal.y)) {
        invalid_fragments++;
      }
      else {
        const double surface_pressure =
          embed_interpolate (point, pressure_field, barycentre) - p_reference;
        if (!isfinite(surface_pressure))
          invalid_fragments++;
        else {
          pressure_x += dA*surface_pressure*normal.x;
          area_sum += dA;
          max_surface_pressure = max(max_surface_pressure, surface_pressure);
        }

        if (constant(viscosity.x) != 0.) {
          double metric_viscosity_sum = 0., face_metric_sum = 0.;
          foreach_dimension() {
            metric_viscosity_sum += viscosity.x[] + viscosity.x[1];
            face_metric_sum += fm.x[] + fm.x[1];
          }

          if (!(face_metric_sum > 0.) || !isfinite(face_metric_sum) ||
              !isfinite(metric_viscosity_sum))
            invalid_fragments++;
          else {
            const double local_viscosity =
              metric_viscosity_sum/face_metric_sum;
            const coord dudn =
              embed_gradient (point, velocity, barycentre, normal);
            const double local_viscous_x =
              -dA*local_viscosity*
              (dudn.x*(sq(normal.x) + 1.) +
               dudn.y*normal.x*normal.y);

            if (!isfinite(local_viscosity) || !isfinite(dudn.x) ||
                !isfinite(dudn.y) || !isfinite(local_viscous_x))
              invalid_fragments++;
            else
              viscous_x += local_viscous_x;
          }
        }
      }
    }
  }

  if (fragments < 0.5) {
    invalid_fragments++;
    max_surface_pressure = NAN;
  }

  AxiEmbedImpactForce result = {
    pressure_x, viscous_x, area_sum, max_surface_pressure,
    (int)fragments, (int)invalid_fragments
  };
  return result;
}

/** Area-weighted pressure in fluid cells adjacent to the radial far field. */
static inline AxiTopPressureDiagnostic
axi_top_adjacent_pressure (scalar pressure_field)
{
  double pressure_area = 0., area_sum = 0.;
  double invalid_faces = 0.;

  foreach_boundary (top, reduction(+:pressure_area) reduction(+:area_sum) reduction(+:invalid_faces)) {
    const double dA = 2.*pi*fm.y[0,1]*Delta;
    if (dA > 0.) {
      if (!isfinite(dA) || !isfinite(pressure_field[]))
        invalid_faces++;
      else {
        pressure_area += pressure_field[]*dA;
        area_sum += dA;
      }
    }
  }

  AxiTopPressureDiagnostic result = {
    area_sum > 0. ? pressure_area/area_sum : NAN,
    area_sum,
    (int)invalid_faces + (area_sum > 0. ? 0 : 1)
  };
  return result;
}

#endif
