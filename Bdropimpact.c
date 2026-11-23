/*
 * PRODUCTION DRIVER: 2D AXISYMMETRIC DROP IMPACT ON AN EMBEDDED SOLID.
 *
 * Student editing contract:
 *   1. Change case inputs only in constants.h.
 *   2. Do not tune physics, contact-angle machinery, geometry, AMR ownership,
 *      or output logic in this driver.
 *   3. Basilisk transports the physical VOF field f. The selected embedded
 *      wall model is locked in constants.h.
 *
 * Coordinate convention under axi.h: x is the axial coordinate and y >= 0
 * is the radial coordinate.  The sessile-cap section defined below becomes
 * the rigid axisymmetric solid after revolution about y = 0.  Basilisk's AXI
 * metric is already included in dv(), cm, fm and face fluxes; never multiply
 * solver terms by y a second time.  The factor 2*pi used in diagnostics
 * restores the omitted azimuthal integral and does not alter the equations.
 */
#include "constants.h"
double drop_time_1file, bubble_time_1file, writefile_time_1file, simulation_time_1file; // Per-output-interval wall-clock accounting; not physical simulation time.
double drop_time_total, bubble_time_total, writefile_time_total, simulation_time_total; // Accumulated wall-clock accounting used only in runtime reports.
clock_t simulation_str_time, simulation_end_time; // CPU-clock markers for the interval timing diagnostic.
/* Component-removal helper is intentionally disabled for the gas-cushion
 * model so no resolved liquid or gas component can be deleted. */
#if 0
static void remove_dropletsv(scalar c, bool droplet, double dropsize);
#endif
struct CFDValues cfdbv; // Single runtime parameter record populated from constants.h by numericalmainvalues().
/* Embedded-wall boundary conditions: zero normal velocity enforces
 * impermeability, zero tangential velocity enforces no slip, and zero normal
 * pressure gradient is the native pressure closure at the stationary wall. */
u.n[embed] = dirichlet(0); // No penetration through the embedded solid.
u.t[embed] = dirichlet(0); // No slip along the embedded solid.
p[embed] = neumann(0); // Homogeneous wall-normal pressure gradient.

/* The external-boundary condition f[left]=0 supplies gas-valued cells outside
 * a planar wall. For an arbitrary EMBED geometry, fully solid cs=0 cells are
 * the corresponding solid-side stencil. They carry no physical fluid volume
 * because cm=0, so this does not clip or alter transported liquid in cut/open
 * cells. sf is reset with f because FILTERED two-phase properties use it. */
static void enforce_embedded_gas_cushion (void)
{
	foreach ()
		if (cs[] <= 0.) {
			f[] = 0.;
			sf[] = 0.;
		}
	boundary ({f, sf});
}

/* This overload is registered after native tension.h, so Basilisk executes it
 * after VOF advection and before native curvature and interfacial forcing. */
#if EMBEDDED_GAS_CUSHION_MODEL
event acceleration (i++)
{
	enforce_embedded_gas_cushion();
}
#endif


/* Right and top are open far-field boundaries.  Pressure is referenced to
 * zero there and the normal velocity uses a zero-gradient outflow closure. */
 u.n[right] = neumann(0);   // Axial outflow: zero normal gradient of normal velocity.
 p[right] = dirichlet(0);   // Axial far-field gauge pressure.
 u.n[top] = neumann(0);     // Radial outflow: zero normal gradient of normal velocity.
 p[top] = dirichlet(0);     // Radial far-field gauge pressure.
// The bottom boundary y = 0 retains Basilisk's AXI symmetry condition; it is the rotation axis, not a physical wall.

int main(int argc, char **argv)
{
	simulation_str_time = clock(); // Start wall-clock accounting before case setup.
	simulation_time_1file = 0.0; // Initialize current output-interval runtime.
	writefile_time_1file = 0.0; // Initialize current output-interval I/O time.
	drop_time_1file = 0.0; // Retained legacy timing channel; it does not modify f.
	bubble_time_1file = 0.0; // Retained legacy timing channel; it does not modify f.
	numericalmainvalues(argv, argc, &cfdbv); // Convert the constants.h nondimensional controls into the runtime record.
	;
	size(cfdbv.domainsize); // Set the square meridional computational-domain width.
	origin(0., 0.); // Place the axial/radial domain origin at (x,y) = (0,0).
	int initialgrid = pow(2, LEVELmin); // Convert the minimum refinement level to the base cells per direction.
	init_grid(initialgrid); // Create the initial 2D tree grid; AXI metrics are supplied by the included solver stack.
	;
	char comm[80];
    sprintf (comm, "mkdir -p intermediate");
    system(comm);
	rho1 = cfdbv.rhoL; // Assign nondimensional liquid density to Basilisk phase 1.
	rho2 = cfdbv.rhoG; // Assign nondimensional gas density to Basilisk phase 2.
	mu1 = cfdbv.muL; // Assign nondimensional liquid viscosity to Basilisk phase 1.
	mu2 = cfdbv.muG; // Assign nondimensional gas viscosity to Basilisk phase 2.
	f.sigma = cfdbv.Sigma; // Assign nondimensional surface tension to the sole physical VOF tracer.
	G.x -= 1.0/sq(cfdbv.Froude);  // Apply nondimensional gravity in the negative axial direction: g* = 1/Fr^2.
    Z.x = 0.0; // Keep the reduced-gravity reference acceleration at zero in the axial direction.
    /* The timestep remains Basilisk-controlled.  Only the pressure/projection
     * residual tolerance is overridden here; the commented controls below
     * are documentation and are not active solver settings. */
    //DT = 1.0e-4;          // Inactive: no fixed timestep is imposed.
    //NITERMIN = 1;         // Inactive: Basilisk's default minimum iteration count is retained.
     //NITERMAX = 100;       // Inactive: Basilisk's default maximum iteration count is retained.
	TOLERANCE = 1e-6;       // Absolute multigrid residual tolerance for the pressure projection.
	run(); // Enter Basilisk's event scheduler and advance the 2D AXI solution.
	return 0; // Return success only after the scheduled end condition is reached.
}

event defaults(i = 0)
{
	interfaces = list_add(NULL, f); // Register only f as the transported physical interface; no auxiliary VOF owns transport.
	interfaces = list_add(interfaces,fb); // Bubble tracer
}

event init(i = 0)
{
	const double sessile_solid_center = cfdbv.rs - cfdbv.xm; // Generating-sphere center; the radial center is on the AXI axis y=0.
	const double solid_top = cfdbv.sessileheight; // Runtime validation in constants.h enforces solid_top = sessile_solid_center + rs.
	const double x0 = solid_top + cfdbv.initialdis + 0.50*cfdbv.diameter; // Falling-drop center: solid top + prescribed surface gap + drop radius.
	double Bubtx0 = solid_top + cfdbv.initialdis + cfdbv.diameter * 0.50;// This is center bubble centre same as drop cneter
  //double Bubtx0 = cfdbv.pooldepth + cfdbv.initialdis + cfdbv.diameter - ((cfdbv.dbdelta * cfdbv.diameter)) - 0.50 * (cfdbv.bubblediameter * cfdbv.diameter);
	/* Build cs/fs on the current tree before any AMR operation.  The former
	 * manual refine(..., all, ...) calls invoked refine_embed_linear() while
	 * the coarse embedded fractions were still uninitialised, causing the
	 * embed-tree.h:280 assertion during event init.  The native adapt event
	 * owns subsequent geometry-first refinement after this initialization. */
	solid (cs, fs, sq(x - sessile_solid_center) + sq(y) - sq(cfdbv.rs)); // Keep only the sphere portion in x>=0 as a rigid sessile-cap solid; positive level-set values denote exterior fluid.
	/* Native EMBED cleanup may change both cs and fs.  Complete it before
	 * constructing the combined AXI+EMBED metrics so cm/fm represent the final
	 * cleaned geometry used by the pressure and timestep operators. */
	fractions_cleanup(cs, fs);
	/* Combine native EMBED fractions with native AXI radial metrics.  cm and
	 * fm are metric-weighted fluid measures, not replacements for cs and fs. */
	cm_update(cm, cs, fs); // Update the AXI+EMBED cell metric from cs and fs.
	fm_update(fm, cs, fs); // Update the AXI+EMBED face metrics from cs and fs.
	restriction({cs, fs, cm, fm}); // Restrict geometry and metric fields consistently over the tree hierarchy.
	boundary({cs, fs, cm, fm}); // Synchronize geometry and metric ghost values, including MPI boundaries.

	/* Initialize the liquid sphere using Basilisk's geometric fraction
	 * operator; this creates conservative PLIC-compatible mixed cells. */
	//fraction(f, sq(0.50*cfdbv.diameter) - (sq(x - x0) + sq(y))); // f=1 in liquid, f=0 in gas, and 0<f<1 at the interface.
	fraction(f, min(sq (0.50 * cfdbv.diameter) - (sq (x - x0) + sq (y)), -(sq (0.50 * (cfdbv.bubblediameter * cfdbv.diameter)) - (sq (x - Bubtx0) + sq (y))))); // this is total fraction okay drop + Bubble
#if EMBEDDED_GAS_CUSHION_MODEL
	enforce_embedded_gas_cushion(); // Initialize the arbitrary solid-side stencil as pure gas.
#endif
    //refine(sq(x - x0) + sq(y) + sq(z) < sq(0.50*cfdbv.diameter + cfdbv.refinegap) && sq(x - x0) + sq(y) + sq(z) > sq(0.50*cfdbv.diameter - cfdbv.refinegap) && level < LEVELmax); // refinement along Dorp 
	foreach () {
		u.x[] = -cfdbv.vel*f[]; // Initialize liquid impact velocity in the negative axial direction; gas starts quiescent.
		u.y[] = 0.0; // Initialize radial velocity to zero in both phases.
		/* two-phase.h owns sf when FILTERED is enabled. Give its first
		 * properties() call the same physical initial fraction as f; later
		 * tracer_advection keeps the existing native filtering policy. */
		sf[] = f[]; // Seed the filtered property fraction from the physical VOF field without changing f.
	}
	boundary ({f, sf, u}); // Synchronize initialized phase and velocity fields before the first solver event.
	/* Create one timing log per MPI rank.  These records report computational
	 * cost only and are not CFD validation data. */
	clock_t timestr, timeend;
	timestr = clock(); // Start timing the initial log creation.
	FILE *fp;
	char name[100], tmp[50];
	sprintf(name, FILENAME_DURATION); // Begin the timing filename with the configured base name.
	sprintf(tmp, "-CPU%02d.plt", pid()); // Add the MPI-rank suffix to avoid concurrent file writes.
	strcat(name, tmp); // Complete the rank-local timing filename.
	fp = fopen(name, "w"); // Create or replace this rank's timing log for the current run.
	fprintf(fp, "Variables = Iteration DeltaTime CriticalTime PhysicalTime LastDuration DropDuration BubbleDuration FileDuration CellNumber TotalLastDuration TotalDropDuration TotalBubbleDuration TotalFileDuration\r\nzone\r\n"); // Write the Tecplot-style column contract.
	fclose(fp); // Close the initialization write immediately.
	timeend = clock(); // Stop timing the initial log creation.
	writefile_time_1file += (double)(timeend - timestr) / CLOCKS_PER_SEC; // Accumulate initialization I/O cost in seconds.
}

typedef struct {
  double liquid, total;
} AxiDissipationRates;

/* Cell-integrated Newtonian viscous power for a stationary no-slip wall.
   x is axial and y radial. dv() contains the radial and fluid-volume metrics. */
static inline AxiDissipationRates axi_embed_cell_dissipation
  (Point point, vector velocity, double liquid_fraction,
   double mu_liquid, double mu_gas)
{
  if (cs[] <= 0.)
    return (AxiDissipationRates){0., 0.};

  coord grad_axial, grad_radial;
  double hoop;
  if (cs[] < 1.) {
    coord wall, fluid;
    coord normal = facet_normal (point, cs, fs);
    double alpha = plane_alpha (cs[], normal);
    plane_area_center (normal, alpha, &wall);
    plane_center (normal, alpha, cs[], &fluid);
    normalize (&normal);

    /* Constant wall velocity makes the tangential derivative zero ON the
       facet. Its first-order Taylor expansion approximates the velocity in
       the fluid part of this cell. Include the native fallback coefficient. */
    coord dudn;
    foreach_dimension() {
      double coefficient = 0.;
      double constant_part = dirichlet_gradient (point, velocity.x, cs, normal,
                                                 wall, 0., &coefficient);
      dudn.x = constant_part + coefficient*velocity.x[];
    }
    grad_axial.x = dudn.x*normal.x;
    grad_axial.y = dudn.x*normal.y;
    grad_radial.x = dudn.y*normal.x;
    grad_radial.y = dudn.y*normal.y;

    double distance = Delta*((fluid.x - wall.x)*normal.x +
                             (fluid.y - wall.y)*normal.y);
    double radius = y + Delta*fluid.y;
    hoop = radius > 0. ? dudn.y*distance/radius : grad_radial.y;
  }
  else {
    scalar axial_velocity = velocity.x, radial_velocity = velocity.y;
    foreach_dimension() {
      grad_axial.x = center_gradient (axial_velocity);
      grad_radial.x = center_gradient (radial_velocity);
    }
    hoop = velocity.y[]/y;
  }

  double strain2 = 2.*(sq(grad_axial.x) + sq(grad_radial.y) + sq(hoop))
                   + sq(grad_axial.y + grad_radial.x);
  double measure = 2.*pi*dv()*strain2;
  double liquid = liquid_fraction*mu_liquid*measure;
  return (AxiDissipationRates){liquid, liquid + (1. - liquid_fraction)*mu_gas*measure};
}

/* Diagnostic accumulators only: these variables never feed momentum, VOF,
 * contact-angle, curvature, AMR, or pressure updates. */
double de, dee, se; 
event energy_budgetDCB (i=0;i++) // Evaluate diagnostic energy terms after every iteration; no physical field is modified.
{    
  static FILE * fp = fopen("energy.txt", "w"); // Open the local energy-diagnostic stream.
  double vd = 0, VD = 0; // Instantaneous liquid-only and two-phase viscous-dissipation estimates.
	double ke1 = 0, ke2 = 0; // Axial and radial kinetic-energy integrals over both phases.
	double ke = 0; // Liquid-only kinetic-energy integral.
	double area = 0; // Reconstructed liquid-gas interfacial area.
	double gpe = 0; // Gravity potential-energy accumulator; its source line below is intentionally inactive.
	double kn = (12./((1 - cfdbv.bubblediameter * cfdbv.bubblediameter * cfdbv.bubblediameter)*pi)); // Legacy nondimensional diagnostic normalization; it does not enter the solver.
	double PreFactor = 2*pi; // Restore the azimuthal 2*pi factor omitted by Basilisk's AXI metric convention.
	foreach (reduction(+:ke1) reduction(+:ke2) reduction(+:ke) reduction(+:area) reduction(+:vd) reduction(+:VD) reduction(+:gpe))
	{
		/* AXI dv() already contains the radial metric.  PreFactor restores only
		 * the azimuthal integration; no additional radial factor is permitted. */
		ke1 +=  PreFactor*0.5*rho1*dv()*f[]*sq(u.x[]) + PreFactor*0.5*rho2*dv()*(1 - f[])*sq(u.x[]); // Two-phase axial kinetic energy.
		ke2 +=  PreFactor*0.5*rho1*dv()*f[]*sq(u.y[]) + PreFactor*0.5*rho2*dv()*(1 - f[])*sq(u.y[]); // Two-phase radial kinetic energy.
		ke +=  PreFactor*0.5*rho1*(sq(u.x[]) + sq(u.y[]))*dv()*f[]; // Liquid-only total kinetic energy.
		//gpe +=  PreFactor*((rho1*f[]) + (rho2*(1 - f[])))*dv()*(1.0/sq(cfdbv.Froude))*x ; // Inactive: GPE remains zero, so TEYG is not an active gravity-energy closure.
		if (f[] > 1e-6 && f[] < 1. - 1e-6)
		{
			coord p;
			coord n = mycs (point, f); // Reconstruct the local PLIC interface normal from f.
			double alpha = plane_alpha (f[], n); // Recover the local PLIC plane constant.
			double s = plane_area_center (n, alpha, &p); // Compute the dimensionless interface segment length in this cell.
			//area = interface_area (f); // Inactive alternative; no global helper result is substituted.
			area +=  PreFactor*s*dv()/Delta; // Convert the segment measure to revolved AXI interfacial area.
		}
		;
		AxiDissipationRates rates = axi_embed_cell_dissipation (point, u, f[], mu1, mu2);
		vd += rates.liquid;
		VD += rates.total;
	}
	/* Normalize instantaneous energies with the existing diagnostic factor. */
	double KE1 = kn*ke1; // Normalized two-phase axial kinetic energy.
	double KE2 = kn*ke2; // Normalized two-phase radial kinetic energy.
	double KE = KE1+KE2; // Normalized total two-phase kinetic energy.
	double KEE = kn*ke; // Normalized liquid-only kinetic energy.
	se = area*cfdbv.Sigma; // Surface energy from reconstructed area and nondimensional surface tension.
	double SE = kn*se; // Normalized surface energy.
	/* Seed the first sample, then integrate the rates over each elapsed interval. */
	static bool have_previous_rates = false;
	static double previous_t, previous_vd, previous_VD;
	if (have_previous_rates) {
		double elapsed = t - previous_t;
		de += 0.5*(previous_vd + vd)*elapsed;
		dee += 0.5*(previous_VD + VD)*elapsed;
	}
	previous_t = t;
	previous_vd = vd;
	previous_VD = VD;
	have_previous_rates = true;
	double DE = kn*de; // Normalized cumulative liquid-only dissipation.
	double DEE = kn*dee; // Normalized cumulative two-phase dissipation.
	double GPE = kn*gpe; // Always zero while the GPE accumulation line remains commented out.
	;
	double TENG = KE+SE+DEE ; // Reported kinetic + surface + cumulative two-phase dissipation diagnostic.
	double TEYG = KE+SE+DEE+GPE; // Numerically equals TENG because GPE accumulation is inactive.
	;
	fprintf (fp, "%d  %f  %.10f  %.10f  %.10f  %.10f  %.10f  %.10f  %.10f  %.10f  %.10f  %.10f  %.10f\r\n", i, t, KE1, KE2, KE, KEE, SE, DE, DEE, GPE, TENG, TEYG, area); // Append one source-derived diagnostic row.
	fflush (fp); // Make each completed row visible without waiting for program termination.
}

/* INACTIVE LEGACY BLOCK: enclosed by this comment and therefore neither
 * compiled nor executed.  It is not part of the production output contract.
event spreading_diameterA (i = 0) 
{
  FILE *fp;
  fp = fopen ("spreading.txt", "w");
  fprintf (fp, "time  Maximum diameter  \r\n");   // First is time , second is dimater okay  
  fclose (fp); 
}

event spreading_diameterB  (t += SAVE_FILE_EVERY) // This is event for calcualting spearding diamter okay dont chnage any thing okay  
{
  FILE *fp;
  fp = fopen ("spreading.txt", "a");
  scalar pos[];
  position (f, pos, {0,1});
  fprintf (fp, "%g %g\n", t, statsf(pos).max);
  fclose (fp);
}*/


event adapt(i++)
{
	double refine[4]; // Absolute wavelet tolerances for the four fields in REFINE_VAR.
	refine[0] = pow(10.0, REFINE_VALUE_0); // Convert the first base-10 tolerance exponent from constants.h.
	refine[1] = pow(10.0, REFINE_VALUE_1); // Convert the second base-10 tolerance exponent from constants.h.
	refine[2] = pow(10.0, REFINE_VALUE_2); // Convert the third base-10 tolerance exponent from constants.h.
	refine[3] = pow(10.0, REFINE_VALUE_3); // Convert the fourth base-10 tolerance exponent from constants.h.
	/* For AXI+EMBED adaptation, the geometry fractions must be refined before
	 * the metric fields.  Preserve every field in Basilisk's normal update
	 * list, but place cs/fs ahead of cm/fm as required by the native tree
	 * refinement contract. */
	scalar * adapt_fields = list_copy({cs, fs, cm, fm}); // Put geometry before metrics in the AMR update list.
	for (scalar s in all)
		adapt_fields = list_add(adapt_fields, s); // Preserve every other registered Basilisk field during adaptation.
	adapt_wavelet(REFINE_VAR, (double[]){refine[0], refine[1], refine[2], refine[3]}, maxlevel = LEVELmax, minlevel = LEVELmin, list = adapt_fields); // Adapt only between the user-controlled minimum and maximum levels.
	free(adapt_fields); // Release the temporary field list after adapt_wavelet returns.
	/* The inherited centered.h adapt handler now performs its native cs/fs
	 * cleanup, closed-face uf cleanup and first properties refresh. */
}

event synchronize_axi_embed_after_adapt(i++, last)
{
	/* The rigid target is fixed and analytic. Reconstruct cs/fs on the final
	 * adapted tree so every MPI rank derives shared embedded faces from the
	 * same level set before rebuilding the combined AXI+EMBED metrics. */
	const double solid_center = cfdbv.rs - cfdbv.xm;
	solid(cs, fs, sq(x - solid_center) + sq(y) - sq(cfdbv.rs));
	fractions_cleanup(cs, fs);
	restriction({cs, fs});
	boundary({cs, fs});
	cm_update(cm, cs, fs);
	fm_update(fm, cs, fs);
	restriction({cm, fm});
	boundary({cm, fm});
#if EMBEDDED_GAS_CUSHION_MODEL
	enforce_embedded_gas_cushion(); // Remove AMR prolongation values only from newly solid cells.
#endif
	/* centered.h closes uf on newly solid faces locally.  Repeat the operation
	 * on the final synchronized geometry and exchange uf across MPI ranks. */
	foreach_face()
		if (uf.x[] && !fs.x[])
			uf.x[] = 0.;
	boundary((scalar *){uf});
	/* two-phase properties depend on the metric-weighted cm/fm fields and must
	 * be refreshed after their final post-AMR reconstruction. */
	event("properties");
}

/* Droplet and bubble removal events are intentionally compiled out. */
#if 0
event drop_remove(i += REMOVE_DROP_PERIOD)
{
	clock_t timestr, timeend;
	timestr = clock();
	if (t > cfdbv.timecontact)
	{
		remove_dropletsv(f, 1, REMOVE_DROP_SIZE);
	}
	timeend = clock();
	drop_time_1file += (double)(timeend - timestr) / CLOCKS_PER_SEC;
}

event drop_bubble(i += REMOVE_BUBBLE_PERIOD)
{
	clock_t timestr, timeend;
	timestr = clock();
	if (t > cfdbv.timecontact)
	{
		remove_dropletsv(f, 0, REMOVE_BUBBLE_SIZE);
	}
	timeend = clock();
	bubble_time_1file += (double)(timeend - timestr) / CLOCKS_PER_SEC;
}
#endif 


event showiteration(i++)
{
	/* Rank zero prints a compact progress identifier.  tb means time before
	 * nominal contact and ta means time after nominal contact. */
	switch (pid())
	{
	case 0:
	{
		char name[500], tmp[100];
		if (t - cfdbv.timecontact < 0.0)
			sprintf(name, "i%05d_dt%.2e_tb%.3f_P%02d", i, dt, t - cfdbv.timecontact, (int)(100.0 * t / MAX_TIME));
		else
		sprintf(name, "i%05d_dt%.2e_ta%.3f_P%02d", i, dt, t - cfdbv.timecontact, (int)(100.0 * t / MAX_TIME));
		sprintf(tmp, "_Re%.5f_We%.5f", (double)cfdbv.Reynolds, (double)cfdbv.Weber);
		//sprintf(tmp, "_Re%d_We%d", (int)cfdbv.Reynolds, (int)cfdbv.Weber);
		strcat(name, tmp);
#if AXI
		sprintf(tmp, "_AXI"); // Active production geometry label for this driver.
		strcat(name, tmp); // Append the compile-time geometry label.
#else
#if dimension == 3
		sprintf(tmp, "_3D"); // Generic inactive fallback; this project does not authorize a 3D run.
		strcat(name, tmp);
#else
		sprintf(tmp, "_2D");
		strcat(name, tmp);
#endif
#endif
		sprintf(tmp, "_L%02d%02d", LEVELmin, LEVELmax);
		strcat(name, tmp);
		printf("%s\r\n", name);
	}
	}
}

event end(t = cfdbv.timecontact + cfdbv.timeend)
{
	/* Write a rank-local completion marker at the configured final physical
	 * time.  Its presence proves event completion, not physical validation. */
	FILE *fp;
	char name[500], tmp[100];
	sprintf(name, FILENAME_ENDOFRUN);
	sprintf(tmp, "-CPU%02d.txt", pid());
	strcat(name, tmp);
	fp = fopen(name, "w");
	fprintf(fp, "SimulationTime %e\r\nRemoveDropTime %e\r\nRemoveBubbleTime %e\r\nWriteFileTime %e\r\n", simulation_time_total, drop_time_total, bubble_time_total, writefile_time_total);
	fclose(fp);
}

event outputfiles (t += SAVE_FILE_EVERY)//remaining the beginning time
{
	/* Output cadence is controlled only by SAVE_FILE_EVERY in constants.h.
	 * This event copies pressure for dumping, writes native Basilisk snapshots,
	 * and updates wall-clock logs; it performs no clipping or field repair. */
	clock_t timestr, timeend;
	timestr = clock();
	static FILE *fp;
	char name[500], tmp[100];
	;
	foreach ()
	{
		if (fb[] < 0.0)
			fb[] = 0.0;
		else if (fb[] > 1.0)
			fb[] = 1.0;
	};
	foreach() {
		pressure[] = p[]; // Copy the solved pressure into the dump-enabled diagnostic scalar.
	}
	p.nodump = false;
	sprintf (name, "intermediate/snapshot-%5.4f", t);
	dump(file = name);
	;
	dump(file = FILENAME_LASTFILE); // Refresh the configured latest-state restart file.
	timeend = clock();
	writefile_time_1file += (double)(timeend - timestr) / CLOCKS_PER_SEC;
	int cellnumber = 0; // Count local leaf cells for the rank-local timing report.
    foreach (noauto) {
        cellnumber++;
    }
	simulation_end_time = clock();
	double estimatetimeleft;
	char LDc[100], TDc[100], ETLc[100];
	simulation_time_1file = (double)(simulation_end_time - simulation_str_time) / CLOCKS_PER_SEC;
	;
	simulation_time_total += simulation_time_1file;
	bubble_time_total += bubble_time_1file;
	drop_time_total += drop_time_1file;
	writefile_time_total += writefile_time_1file;
	if (t == 0.0)
		estimatetimeleft = 0.0; // ETA is undefined at t=0; report zero instead of dividing by zero.
	else
		estimatetimeleft = simulation_time_total * (cfdbv.timecontact + cfdbv.timeend) / t - simulation_time_total; // Linear wall-clock extrapolation; not a solver timestep estimate.
	timecalculation(simulation_time_1file, LDc);
	timecalculation(simulation_time_total, TDc);
	timecalculation(estimatetimeleft, ETLc);
	sprintf(name, FILENAME_DURATION);
	sprintf(tmp, "-CPU%02d.plt", pid());
	strcat(name, tmp);
	fp = fopen(name, "a");
	//fprintf(fp, "Variables = Iteration DeltaTime CriticalTime PhysicalTime LastDuration DropDuration BubbleDuration FileDuration CellNumber TotalLastDuration TotalDropDuration TotalBubbleDuration TotalFileDuration\r\nzone\r\n");
	fprintf(fp, "%d %e %e %e %e %e %e %e %d %e %e %e %e\r\n", i, dt, t - cfdbv.timecontact, t, simulation_time_1file, drop_time_1file, bubble_time_1file, writefile_time_1file, cellnumber, simulation_time_total, drop_time_total, bubble_time_total, writefile_time_total);
	fclose(fp);
	simulation_str_time = clock();
	simulation_time_1file = 0.0;
	writefile_time_1file = 0.0;
	drop_time_1file = 0.0;
	bubble_time_1file = 0.0;
	;
	switch (pid())
	{
	case 0:
	{
		printf("\r\nData Files are Written!\r\nDuration Last: %s\r\nTotal: %s, Time Left: %s\r\n\r\n", LDc, TDc, ETLc);
		break;
	}
	}
}


/* Historical component-removal implementation, retained but disabled. */
#if 0
static void remove_dropletsv(scalar c, bool droplet, double dropsize)
{
	int j, n;
	scalar m[];
	const double THR = R_VOFLIMIT; //THRESHOLD
	const double delta = cfdbv.domainsize / pow(2.0, LEVELmax);
	double realD2;
	;
	foreach (noauto)
		m[] = (droplet ? (c[] > THR) : (c[] < (1. - THR)));
	n = tag(m);
	double v[n];
	for (j = 0; j < n; j++)
		v[j] = 0.0;
	foreach_leaf()
	{
		if (m[] > 0)
		{
			j = m[] - 1;
#if dimension == 3
			v[j] += Delta * Delta * Delta * (droplet ? c[] : 1. - c[]); // c[]
#else
			v[j] += Delta * Delta * (droplet ? c[] : 1. - c[]); // c[]
#endif
			//v[j] += dv()*c[];
		}
	}
#if _MPI
	MPI_Allreduce(MPI_IN_PLACE, v, n, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#endif
	foreach (noauto)
	{
		if (m[] > 0)
		{
			j = m[] - 1;
			realD2 = 4.0 * v[j] / R_PI;
			if (realD2 < dropsize * dropsize * delta * delta) //// we compute the equivalent radius of the droplet. http://basilisk.fr/src/test/spurious.c 
				c[] = droplet ? 0. : 1.;
		}
	}
	boundary({c});
}
#endif
