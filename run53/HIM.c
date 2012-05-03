/**********************************************************************
 *                   GNU General Public License                       *
 * This file is a part of HIM.                                        *
 *                                                                    *
 * HIM is free software; you can redistribute it and/or modify it and *
 * are expected to follow the terms of the GNU General Public License *
 * as published by the Free Software Foundation; either version 2 of  *
 * the License, or (at your option) any later version.                *
 *                                                                    *
 * HIM is distributed in the hope that it will be useful, but WITHOUT *
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY *
 * or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public   *
 * License for more details.                                          *
 *                                                                    *
 * For the full text of the GNU General Public License,               *
 * write to: Free Software Foundation, Inc.,                          *
 *           675 Mass Ave, Cambridge, MA 02139, USA.                  *
 * or see:   http://www.gnu.org/licenses/gpl.html                     *
 **********************************************************************/

/********+*********+*********+*********+*********+*********+*********+*
 *                                                                    *
 *                   The Hallberg Isopycnal Model                     *
 *                               HIM                                  *
 *                                                                    *
 *  By Robert Hallberg                        rwh@gfdl.noaa.gov       *
 *                                                                    *
 *  This file worked on June 1992 - June 2002.                        *
 *                                                                    *
 *    This program (HIM) simulates the ocean by numerically solving   *
 *  the Boussinesq primitive equations in isopycnal vertical coord-   *
 *  inates and general orthogonal horizontal coordinates.  These      *
 *  equations are horizontally discretized on an Arakawa C-grid.      *
 *  There are a range of options for the physical parameterizations,  *
 *  from those most appropriate to highly idealized models for studies*
 *  of geophysical fluid dynamics to a rich suite of processes        *
 *  appropriate for realistic ocean simulations.  The thermodynamic   *
 *  options range from an adiabatic model with fixed density layers   *
 *  to a model with temperature and salinity as state variables and   *
 *  using a full nonlinear equation of state.  The uppermost few      *
 *  layers may be used to describe a bulk mixed layer, including the  *
 *  effects of penetrating shortwave radiation.  Either a split-      *
 *  explicit time stepping scheme or a non-split scheme may be used   *
 *  for the dynamics, while the time stepping may be split (and use   *
 *  different numbers of steps to cover the same interval) for the    *
 *  forcing, the thermodynamics, and for the dynamics.  Most of the   *
 *  numerics are second order accurate in space.  HIM can run with an *
 *  absurdly thin minimum layer thickness.                            *
 *                                                                    *
 *    Details of the numerics and physical parameterizations are      *
 *  provided in the appropriate source files.  Most of the available  *
 *  options are selected by the settings in init.h, although some     *
 *  (such as the equation of state) are selected by specifying which  *
 *  file to use in the make file (Makefile).                          *
 *                                                                    *
 *    The present file (HIM.c) contains the main time stepping loops. *
 *  One time integration option for the dynamics uses a split explicit*
 *  time stepping scheme to rapidly step the barotropic pressure and  *
 *  velocity fields. The barotropic velocities are averaged over the  *
 *  baroclinic time step before they are used to advect thickness     *
 *  and determine the baroclinic accelerations. At the end of every   *
 *  time step, the free surface height perturbation is determining    *
 *  by adding up the layer thicknesses; this perturbation is used to  *
 *  drive the free surface heights from the barotropic calculation    *
 *  and from the sum of the layer thicknesses toward each other over  *
 *  subsequent time steps. The barotropic and baroclinic velocities   *
 *  are synchronized as part of the vertical viscosity algorithm and  *
 *  be recalculating. the barotropic velocities from the baroclinic   *
 *  velocities each time step. This scheme is described in Hallberg,  *
 *  1997, J. Comp. Phys. 135, 54-65.                                  *
 *                                                                    *
 *    The other time integration option uses a non-split time stepping*
 *  scheme based on the 3-step third order Runge-Kutta scheme         *
 *  described in Matsuno, 1966, J. Met. Soc. Japan,  44, 85-88. *                                                                    *
 *    There are a range of closure options available.  Horizontal     *
 *  velocities are subject to a combination of horizontal biharmonic  *
 *  and Laplacian friction (based on a stress tensor formalism) and a *
 *  vertical Fickian viscosity (perhaps using the kinematic viscosity *
 *  of water).  The horizontal viscosities may be constant, spatially *
 *  varying or may be dynamically calculated using Smagorinsky's      *
 *  approach.  A diapycnal diffusion of density and thermodynamic     *
 *  quantities is also allowed, but not required, as is horizontal    *
 *  diffusion of interface heights (akin to the Gent-McWilliams       *
 *  closure of geopotential coordinate models).  The diapycnal mixing *
 *  may use a fixed diffusivity or it may use the shear Richardson    *
 *  number dependent closure described in Hallberg (MWR, 2000).       *
 *  When there is diapycnal diffusion, it applies to momentum as well.*
 *  As this is in addition to the vertical viscosity, the vertical    *
 *  Prandtl always exceeds 1.                                         *
 *                                                                    *
 *    HIM has a number of noteworthy debugging capabilities.          *
 *  Excessively large velocities are truncated and HIM will stop      *
 *  itself after a number of such instances to keep the model from    *
 *  crashing altogether, and the model state is output with a         *
 *  reported time of 9.9e9.  This is useful in diagnosing failures,   *
 *  or (by accepting some truncations) it may be useful for getting   *
 *  the model past the adjustment from an ill-balanced initial        *
 *  condition.  In addition, all of the accelerations in the columns  *
 *  with excessively large velocities may be directed to a text file. *
 *  Parallelization errors may be diagnosed with the CHECK_PARALLEL   *
 *  option, whereby ostensibly identical model incarnations are run   *
 *  simultaneously on one and multiple processors and any differences *
 *  are reported.                                                     *
 *                                                                    *
 *    About 35 other files of source code and 6 header files must     *
 *  be used to make the model work.  A makefile is included to        *
 *  make compiling easy.  Type "make HIM" to compile.  Some run       *
 *  time input is required, but this is prompted for.  It may be      *
 *  convenient to direct a small file containing the needed inform-   *
 *  ation to standard input, and direct the output to another file.   *
 *                                                                    *
 *    The ~35 source files contain the following subroutines:         *
 *  HIM.c:                                                            *
 *    step_HIM steps HIM over a specified interval of time.           *
 *    HIM_initialize calls initialize and does other initialization   *
 *      that does not warrant user modification.                      *
 *    set_restart_fields is used to specify those fields that are     *
 *      written to and read from the restart file.                    *
 *    calculate_surface_state determines the surface (mixed layer)    *
 *      properties of the current model state and packages pointers   *
 *      to these fields into an exported structure.                   *
 *  HIM_checkparallel.c:                                              *
 *   (HIM_checkparallel.c is identical to HIM.c, except that it has   *
 *    a number of diagnostic calls to assess whether a parallel run   *
 *    is giving identical results to a single processor run.)         *
 *  HIM_driver.c:                                                     *
 *    main is where HIM starts.  Inside of main are the calls that    *
 *      set up the run, step the model, and orchestrate output and    *
 *      normal termination of the run.                                *
 *                                                                    *
 *                                                                    *
 *     THE FOLLOWING FILES ARE WHERE INITIAL CONDITIONS, FORCING,     *
 *     AND DOMAIN PROPERTIES ARE PRINCIPALLY SPECIFIED.               *
 *                                                                    *
 *  initialize.c:                                                     *
 *    initialize does just that to all of the fields that are needed  *
 *      to specify the initial conditions of the model.  initialize   *
 *      calls a number of other subroutines in initialize.c, each of  *
 *      which initializes a single field (or a few closely related    *
 *      fields) that are indicated by the subroutine name.            *
 *  initialize_output.c:                                              *
 *    USER_set_output specifies what fields are output into which     *
 *      files and when, and whether they are averaged in time.        *
 *    USER_set_diagnostics specifies what diagnostic quantities are   *
 *      calculated.                                                   *
 *    GetInputLines gets 4 controlling inputs from stdin, including   *
 *      the directories for I/O and the parameter specification file. *
 *    write_grid_file writes out a file describing the model grid.    *
 *  surface_forcing.c:                                                *
 *    set_forcing sets the current values of surface forcing fields.  *
 *    wind_forcing sets the current surface wind stresses.            *
 *    buoyancy_forcing sets the current surface heat, fresh water,    *
 *      buoyancy or other appropriate tracer fluxes.                  *
 *    set_forcing_output sets up the output of any forcing fields.    *
 *    average_forcing accumulates time averages of indicated forcing  *
 *      fields.                                                       *
 *    register_forcing_restarts is used to specify the forcing-related*
 *      fields that are written to and read from the restart file.    *
 *  set_metrics.c:                                                    *
 *    set_metrics calculates the horizontal grid spacings and related *
 *      metric fields, along with the gridpoint locations.            *
 *    initialize_masks initializes the land masks.                    *
 *                                                                    *
 *     THE FOLLOWING FILES CONTAIN THE PRINCIPAL DYNAMIC ROUTINES.    *
 *                                                                    *
 *  CoriolisAdv.c:                                                    *
 *    CorAdCalc calculates the Coriolis and advective accelerations.  *
 *  PressureForce.c:                                                  *
 *    PressureForce calculates the pressure acceleration.             *
 *    register_compress is used to specify the reference profile of   *
 *      potential temperature and salinity that is used to compensate *
 *      for compressibility.                                          *
 *    uncompress_e_rho makes internally consistent changes to profiles*
 *      of interface height and density to offset compressibility and *
 *      to minimize the non-solenoidal pressure gradient term.        *
 *  continuity.c, continuity_mpdata.c, or continuity_FCT.c:           *
 *    continuity time steps the layer thicknesses.                    *
 *  barotropic.c:                                                     *
 *    btstep time steps the linearized barotropic equations for use   *
 *      with the split explicit time stepping scheme.                 *
 *    btcalc calculates the barotropic velocities from the layer      *
 *      velocities.                                                   *
 *    barotropic_init initializes several split-related variables and *
 *      calculates several static quantities for use by btstep.       *
 *    register_barotropic_restarts indicates those time splitting-    *
 *      related fields that are to be in the restart file.            *
 *  hor_visc.c:                                                       *
 *    horizontal_viscosity calculates the convergence of momentum     *
 *      due to Laplacian or biharmonic horizontal viscosity.          *
 *    set_up_hor_visc calculates combinations of metric coefficients  *
 *      and other static quantities used in horizontal_viscosity.     *
 *  vertvisc.c:                                                       *
 *    vertvisc changes the velocity due to vertical viscosity,        *
 *      including application of a surface stress and bottom drag.    *
 *    set_viscous_BBL determines the bottom boundary layer thickness  *
 *      and viscosity according to a linear or quadratic drag law.    *
 *  thickness_diffuse.c:                                              *
 *    thickness_diffuse moves fluid adiabatically to horizontally     *
 *      diffuse interface heights.                                    *
 *                                                                    *
 *                                                                    *
 *     THE FOLLOWING FILES CONTAIN THE THERMODYNAMIC ROUTINES.        *
 *                                                                    *
 *  diabatic_driver.c:                                                *
 *    diabatic orchestrates the calculation of vertical advection and *
 *      diffusion of momentum and tracers due to diapycnal mixing and *
 *      mixed layer (or other diabatic) processes.  mixedlayer,       *
 *      Calculate_Entrainment, apply_sponge, and any user-specified   *
 *      tracer column physics routines are all called by diabatic.    *
 *  diabatic_entrain.c:                                               *
 *    Calculate_Entrainment calculates the diapycnal mass fluxes due  *
 *      to interior diapycnal mixing processes, which may include a   *
 *      Richardson number dependent entrainment.                      *
 *    Calculate_Rino_flux estimates the Richardson number dependent   *
 *      entrainment in the absence of interactions between layers,    *
 *      from which the full interacting entrainment can be found.     *
 *    Estimate_u_h estimates what the velocities at thickness points  *
 *      will be after entrainment.                                    *
 *  mixed_layer.c:                                                    *
 *    mixed_layer implements a bulk mixed layer, including entrainment*
 *      and detrainment, related advection of dynamically active      *
 *      tracers, and buffer layer splitting.  The bulk mixed layer    *
 *      may consist of several layers.                                *
 *  tracer.c:                                                         *
 *    register_tracer is called to indicate a field that is to be     *
 *      advected by advect_tracer and diffused by tracer_hordiff      *
 *    advect_tracer does along-isopycnal advection of tracer fields.  *
 *    tracer_hordiff diffuses tracers along isopycnals.               *
 *    register_tracer_init_fn registers a user-specified tracer       *
 *      initialization subroutine.                                    *
 *    call_tracer_init_fns calls any user-specified tracer initial-   *
 *      ization subroutines that have been registered.                *
 *    register_tracer_column_fn registers a user-specified tracer     *
 *      column processes subroutine.                                  *
 *    call_tracer_column_fns calls any user-specified tracer column   *
 *      processes subroutines that have been registered.              *
 *  sponge.c:                                                         *
 *    apply_sponge damps fields back to reference profiles.           *
 *    initialize_sponge stores the damping rates and allocates the    *
 *      memory for the reference profiles.                            *
 *    set_up_sponge_field registers reference profiles and associates *
 *      them with the fields to be damped.                            *
 *  eqn_of_state.c, eqn_of_state_linear.c, or eqn_of_state_UNESCO.c:  *
 *    calculate_density calculates a list of densities at given       *
 *      potential temperatures, salinities and pressures.             *
 *    calculate_density_derivs calculates a list of the partial       *
 *      derivatives with temperature and salinity at the given        *
 *      potential temperatures, salinities and pressures.             *
 *    calculate_compress calculates a list of the compressibilities   *
 *      (partial derivatives of density with pressure) at the given   *
 *      potential temperatures, salinities and pressures.             *
 *    calculate_2_densities calculates a list of the densities at two *
 *      specified reference pressures at the given potential          *
 *      temperatures and salinities.                                  *
 *  fit_compressibility.c:                                            *
 *    fit_compressibility determines the best fit of compressibility  *
 *      with pressure using a fixed 5-coefficient functional form,    *
 *      based on a provided reference profile of potential temperature*
 *      and salinity with depth.  This fit is used in PressureForce.  *
 *                                                                    *
 *                                                                    *
 *     THE FOLLOWING FILES CONTAIN INFRASTRUCTURAL ROUTINES.          *
 *                                                                    *
 *  HIM_restart.c:                                                    *
 *    save_restart saves a restart file (or multiple files if they    *
 *      would otherwise be too large).                                *
 *    register_restart_field is called to specify a field that is to  *
 *      written to and read from the restart file.                    *
 *    restore_state reads the model state from restart or other files.*
 *    query_initialized indicates whether a specific field or all     *
 *      restart fields have been read from the restart files.         *
 *  HIM_parser.c:                                                     *
 *    HIM_parser parses a parameter specification file with the same  *
 *      format as init.h to (possibly) change parameters at run time. *
 *  HIM_sum_output.c:                                                 *
 *    write_energy writes the layer energies and masses and other     *
 *      spatially integrated quantities and monitors CPU time use.    *
 *    depth_list_setup generates a list of the volumes of fluid below *
 *      various depths.                                               *
 *  HIM_field_output.c:                                               *
 *    specify_saves_file sets up a file for output, including when and*
 *      how frequently the file is to be written.                     *
 *    set_up_save_field sets up a field to be saved and also perhaps  *
 *      to be averaged in time.                                       *
 *    save_fields saves all fields that are to be saved within the    *
 *      given time interval and indicates the next scheduled time     *
 *      that a field is to be saved.                                  *
 *    enable_averaging enables averaging for a time interval.         *
 *    disable_averaging disables the accumulation of averages.        *
 *    query_averaging_enabled indicates whether averaging is          *
 *      currently enabled.                                            *
 *    get_field_ref returns a reference integer for a field that is   *
 *      to be averaged, or 0 if it is not to be averaged.             *
 *    average_field accumulates the average of a field.               *
 *  timetype.c:   (Partial C implementation of the FMS time manager.) *
 *    t_ge, t_gt, t_le, t_lt, t_eq, and t_ne make logical comparisons *
 *      between times.                                                *
 *    t_plus adds two times together.                                 *
 *    t_minus finds the absolute difference between times. (|t1-t2|)  *
 *    t_mult multiplies a time by an integer.                         *
 *    t_div returns the real ratio of two times.                      *
 *    t_div_int divides a time by an integer.                         *
 *    t_interval finds the signed, real ratio of the interval between *
 *      two times to a third time.  (i.e. (t1-t2)/t3)                 *
 *    double_to_time converts a real number to a time.                *
 *    time_to_double converts a time to a real number.                *
 *  HIM_io.c:         (Input/Output utility subroutines.)             *
 *    create_file creates a new file, set up structures that are      *
 *      needed for subsequent output, and write the coordinates.      *
 *    reopen_file reopens an existing file for writing and set up     *
 *      structures that are needed for subsequent output.             *
 *    open_input_file opens the indicated file for reading only.      *
 *    close_file closes an open file.                                 *
 *    synch_file flushes the buffers for a file, completing all       *
 *      pending output.                                               *
 *    write_field writes a field to an open file.                     *
 *    write_time writes a value of the time axis to an open file.     *
 *    Read_Field provides a simpler interface to read_field.          *
 *    read_field reads a field from an open file.                     *
 *    read_time reads a time from an open file.                       *
 *    name_output_file provides a name for an output file based on a  *
 *      name root and the time of the output.                         *
 *    find_input_file finds a file that has been previously written by*
 *      HIM and named by name_output_file and opens it for reading.   *
 *    handle_error writes an error code and quits.                    *
 *  HIM_parallel.c:    (Parallelization utility subroutines.)         *
 *    set_up_parallel starts parallel runs and calculates subdomain   *
 *      sizes and neighboring domains.                                *
 *    pass_var passes a 3-D variable to neighboring processors or     *
 *      applies corresponding boundary conditions.                    *
 *    pass_var_2d passes a 2-D variable to neighboring processors or  *
 *      applies corresponding boundary conditions.                    *
 *    sum_fields sums fields across all processors.                   *
 *    sum_int_fields sums integer fields across all processors.       *
 *    collect_double_fields collects a full layer of a field to       *
 *      processor 0.                                                  *
 *    spread_string distributes a character string from processor 0.  *
 *    spread_double_vector distributes copies of a 1-D array from     *
 *      processor 0.                                                  *
 *    check_field checks whether a field is identical on different    *
 *      processors, and reports any differences.  This is useful for  *
 *      checking parallelization or for debugging coding changes that *
 *      are not expected to change answers.                           *
 *    quit causes all processors to stop execution.                   *
 *                                                                    *
 *                                                                    *
 *     THE FOLLOWING FILES CONTAIN PURELY DIAGNOSTIC ROUTINES.        *
 *                                                                    *
 *  diagnostics.c:                                                    *
 *    calculate_diagnostic_fields is used to calculate several        *
 *      diagnostic fields that are not naturally calculated elsewhere.*
 *    register_time_deriv is used to register the information needed  *
 *      for diagnostically calculating a time derivative.             *
 *    calculate_derivs calculates any registered time derivatives.    *
 *  PointAccel.c:                                                     *
 *    write_u_accel writes a long list of zonal accelerations and     *
 *      related quantities for one column out to a file.  This is     *
 *      typically called for diagnostic purposes from vertvisc when a *
 *      zonal velocity exceeds the specified threshold.               *
 *    write_v_accel writes a long list of meridional accelerations and*
 *      related quantities for one column out to a file.  This is     *
 *      typically called for diagnostic purposes from vertvisc when a *
 *      meridional velocity exceeds the specified threshold.          *

 *                                                                    *
 *    In addition there are 5 header files:                           *
 *  init.h sets various constants and parameters for the simulation.  *
 *  HIM.h contains the subroutine prototypes and the definitions of   *
 *    the following 6 structures:                                     *
 *      params (HIM_params) contains changeable simulation parameters.*
 *      forcing (fluxes) contains pointers to surface forcing fields. *
 *      surface (state) contains pointers to ocean surface fields.    *
 *      diag_fld (diag) contains pointers to diagnostic quantities.   *
 *      thermo_var_ptrs (tv) contains pointers to thermodynamic       *
 *        fields, such as potential temperature and salinity.         *
 *      bt_vars_ptrs (bt) contains pointers to variables related to   *
 *        the baroclinic/barotropic split time stepping.              *
 *  HIM_io.h contains the subroutine prototypes and structure         *
 *    definitions for I/O related subroutines.                        *
 *  metrics.h contains the descriptions for a number of metric terms. *
 *  hor_visc.h contains the descriptions for a number of metric-      *
 *    related fields that are only used by the horizontal viscosity.  *
 *                                                                    *
 *                                                                    *
 *    Most simulations can be set up by modifying only the files      *
 *  init.h, initialize.c, and surface_forcing.c.  In addition, the    *
 *  file initialize_output.c will commonly be modified to tailor the  *
 *  output to the needs of the question at hand.  These altered files *
 *  and the make file might reside in the same directory as the       *
 *  executable file.  All of the other (unaltered) source code should *
 *  probably remain in some central directory.  The makefile will, of *
 *  course, need to reflect the locations of the desired source files.*
 *                                                                    *
 *  Variables written all in capital letters are defined in init.h    *
 *                                                                    *
 *     A small fragment of the grid is shown below:                   *
 *                                                                    *
 *    j+1  x ^ x ^ x   At x:  q, f                                    *
 *    j+1  > o > o >   At ^:  v, PFv, CAv, vh, diffv, tauy, vbt, vhtr *
 *    j    x ^ x ^ x   At >:  u, PFu, CAu, uh, diffu, taux, ubt, uhtr *
 *    j    > o > o >   At o:  h, D, eta, T, S, tr, actflux            *
 *    j-1  x ^ x ^ x                                                  *
 *        i-1  i  i+1                                                 *
 *           i  i+1                                                   *
 *                                                                    *
 *  The boundaries always run through q grid points (x).              *
 *                                                                    *
 ********+*********+*********+*********+*********+*********+*********+*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <init.h>
#define MAINFILE
#include <metrics.h>
#include <HIM.h>
#include <HIM_io.h>
#include <timetype.h>

double u[3][NZ][NYMEM][NXMEM];  /* Zonal velocity, in m s-1.          */
double v[3][NZ][NYMEM][NXMEM];  /* Meridional velocity, in m s-1.     */
double h[3][NZ][NYMEM][NXMEM];  /* Layer thickness, in m.             */
double D[NYMEM][NXMEM];         /* Basin depth, in m.                 */
double CAu[NZ][NYMEM][NXMEM];   /* CAu = f*v - u.grad(u), m s-2.      */
double CAv[NZ][NYMEM][NXMEM];   /* CAv = -f*u - u.grad(v), m s-2.     */
double PFu[NZ][NYMEM][NXMEM];   /* PFu = -dM/dx, in m s-2.            */
double PFv[NZ][NYMEM][NXMEM];   /* PFv = -dM/dy, in m s-2.            */
double uh[NZ][NYMEM][NXMEM];    /* uh = u * h * dy at u grid points.  */
double vh[NZ][NYMEM][NXMEM];    /* vh = v * h * dx at v grid points.  */
                                /* uh and vh are in m3 s-1.           */
double diffu[NZ][NYMEM][NXMEM]; /* Zonal or meridional accelerations  */
double diffv[NZ][NYMEM][NXMEM]; /* due to convergence of the along-   */
                                /* isopycnal stress tensor, in m s-2. */

double uhtr[NZ][NYMEM][NXMEM];/* Zonal & meridional thickness fluxes  */
double vhtr[NZ][NYMEM][NXMEM];/* used to advect tracers, m3 s-1.      */

double hmask[NYMEM][NXMEM];     /* hmask is 0 for land points and 1   */
                                /* for ocean points.                  */
double umask[NYMEM][NXMEM];     /* _mask are 0 for velocity points    */
double vmask[NYMEM][NXMEM];     /* where a point on either side is    */
                                /* land and 1 otherwise.              */
double qmask[NYMEM][NXMEM];     /* qmask is 1 for interior points and */
                                /* 0 for land and boundary points.    */

extern struct params HIM_params;/* HIM_params is a structure that     */
                                /* contains a number of parameters    */
                                /* of the simulation which can be     */
                                /* modified from a parameter          */
                                /* specification file at run time.    */

extern struct diag_fld diag;    /* A structure containing pointers to */
                                /* the diagnostic fields which might  */
                                /* be calculated.                     */

#ifdef SPLIT                                           /* start SPLIT */
double pbce[NZ][NYMEM][NXMEM];/* pbce times s gives the baroclinic    */
                              /* pressure anomaly in each layer due   */
                              /* to free surface height anomalies.    */
                              /* In m s-2.                            */

double u_accel_bt[NYMEM][NXMEM]; /* u_accel_bt and v_accel_bt are the */
double v_accel_bt[NYMEM][NXMEM]; /* difference between the acceler-   */
                              /* ations from the barotropic calcul-   */
                              /* ation and the vertical averages of   */
                              /* the layer accelerations, in m s-2.   */
double ubt_force[NYMEM][NXMEM]; /* ubt_force and vbt_force are the    */
double vbt_force[NYMEM][NXMEM]; /* momentum forcing terms which drive */
                              /* the barotropic model, and these in-  */
                              /* clude both forcing (e.g. wind stress)*/
                              /* and the vertical average of the      */
                              /* accelerations which are not ex-      */
                              /* plictly included in the barotropic   */
                              /* equations.  Both are in m s-2.       */
#endif                                                   /* end SPLIT */

struct thermo_var_ptrs thermovar; /* A structure containing pointers  */
                              /* to an assortment of thermodynamic    */
                              /* fields that may be available, inclu- */
                              /* ding potential temperature, salinity */
                              /* and mixed layer density.             */

struct bt_vars_ptrs bt;      /* A structure containing pointers to a  */
                             /* collection of variables related to the*/
                             /* barotropic-baroclinic split.  Each of */
                             /* these pointers is NULL if the split   */
                             /* time stepping scheme is not used.     */

#ifdef BULKMIXEDLAYER                         /* start BULKMIXEDLAYER */
double Rml[NML+1][NYMEM][NXMEM];/* The mixed and buffer layer potent- */
                              /* ial densities in kg m-3.             */
#endif                                          /* end BULKMIXEDLAYER */

#ifdef TEMPERATURE                               /* start TEMPERATURE */
double T[NZ][NYMEM][NXMEM];   /* Potential temperature, in C.         */
double S[NZ][NYMEM][NXMEM];   /* Salinity in PSU.                     */
# ifdef FRAZIL
double frazil[NYMEM][NXMEM];  /* The accumulated energy required to   */
                              /* keep the temperature above freezing  */
                              /* since the start of the call to       */
                              /* step_HIM.                            */
# endif
#endif                                             /* end TEMPERATURE */

double f[NYMEM][NXMEM];       /* The Coriolis parameter of lat., s-1. */
double g[NZ];                 /* g is the reduced gravity at the      */
                              /* interfaces, in units of m s-2.       */
double Rlay[NZ+1];            /* Rlay is the target potential density */
                              /* (the coordinate variable) of a layer */
                              /* in m s-2.                            */

long ntrunc = 0;              /* ntrunc is the number of times the    */
                              /* velocity has been truncated since    */
                              /* the last call to write_energy.       */

int pe_here = 0;              /* The current processor's label.       */
int X1off = 0, Y1off = 0;     /* The offset in the global indices of  */
                              /* X1 and Y1 on the current processor,  */
                              /* relative to X1 and Y1 on PE 0.       */
int nx = NXTOT+X1-1;          /* The largest index value in the x-    */
int ny = NYTOT+Y1-1;          /* and y- directions of the local       */
                              /* computational domain.                */

double rel_time = 0.0;        /* Relative time in s since the start   */
                              /* of the current execution.            */

void calculate_surface_state(double u[NZ][NYMEM][NXMEM], double v[NZ][NYMEM][NXMEM],
                             double h[NZ][NYMEM][NXMEM], struct surface *state);
static void set_restart_fields(void);

int step_HIM(double time_int, struct forcing fluxes, struct surface *state) {
/* ================================================================== */
/*   This subroutine time steps HIM.                                  */
/* ================================================================== */
/* Arguments: time_int - The interval of time over which to integrate */
/*                       in s.                                        */
/*  (in)      fluxes - A structure containing pointers to any possible*/
/*                     forcing fields.  Unused fields have NULL ptrs. */
/*  (out)     state - A structure containing fields that describe the */
/*                    surface state of the ocean.                     */

  double dt;                  /* The baroclinic time step in s.       */
  double dtnt = 0.0;          /* The elapsed time since updating the  */
                              /* tracers and applying diabatic        */
                              /* processes, in s.                     */
  double dt_pred;             /* The time step for the predictor part */
                              /* of the baroclinic time stepping.     */
  double I_rho0 = 1.0/RHO_0;  /* The inverse of the mean density, in  */
                              /* units of m3 kg-1.                    */
  int m = 0;                  /* The current time level (0, 1, or 2). */
  int mp;                     /* The previous value of m.             */
  int ntstep;                 /* The number of time steps between     */
                              /* tracer updates or diabatic forcing.  */
  int n_max;                  /* The number of steps to take in this  */
                              /* call.                                */
  static int nt=1;            /* The running number of iterations.    */
  int calc_bbl;               /* When 1, the BBL viscosity and thick- */
                              /* ness are calculated within vertvisc. */
                              /* This only applies when BOTTOMDRAGLAW */
                              /* defined.                             */
  int i, j, k, n;

/*   First determine the time step that is consistent with this call. */
/* It is anticipated that the time step will almost always coincide   */
/* with HIM_params.dt.  In addition, ntstep is determined, subject to */
/* the constraint that ntstep cannot exceed n_max.                    */
  n_max = (time_int <= HIM_params.dt) ? 1 :
          (int) ceil(time_int/HIM_params.dt - 0.001);
  dt = time_int / (double) n_max;
#ifdef SPLIT
  dt_pred = dt * HIM_params.be;
#else
  dt_pred = dt / 3.0;
#endif
  ntstep = (HIM_params.ntstep <= 1) ? 1 :
           ((HIM_params.ntstep >= n_max) ? n_max : HIM_params.ntstep);
  calc_bbl = 1;

  if (fluxes.ustar != NULL)
    pass_var_2d(*fluxes.ustar, 1, 1, 1, 0.0, 1, 1, TO_SOUTH+TO_WEST,0);

  if (thermovar.frazil != NULL)
    for (j=Y1;j<=ny;j++) for (i=X1;i<=nx;i++) thermovar.frazil[j][i] = 0.0;


  for (n=1;n<=n_max;n++,nt++) {
    rel_time += dt;

#ifdef SPLIT /* ----------------------------------------- start SPLIT */
/* This section uses a predictor corrector scheme, that is somewhere  */
/* (determined by HIM_params.be) between the forward-backward (be=0.5)*/
/* scheme and the backward Euler scheme (be=1.0) to time step the     */
/* dynamic equations.                                                 */
    m = (int) (nt % 2);
    mp = 1 - m;

    disable_averaging();

    btcalc(u[mp],v[mp],h[mp],mp,1);
    CorAdCalc(u[2],v[2],h[2],uh,vh,CAu,CAv);

# ifndef BEGW
#  define BEGW 0.0
# endif
    if (BEGW == 0.0) enable_averaging(dt,rel_time);
    PressureForce(h[mp],thermovar,bt,PFu,PFv,pbce,1,1);
    disable_averaging();

/* Calculate the momentum forcing terms for the barotropic equations. */
    for (j=Y1;j<=ny;j++) {
      for (i=X1;i<=nx;i++) {
        ubt_force[j][i] = (*(fluxes.taux))[j][i] * I_rho0 * bt.IDatu[j][i] +
          bt.frhatu[0][j][i] * (CAu[0][j][i] + diffu[0][j][i] + PFu[0][j][i]);
        vbt_force[j][i] = (*(fluxes.tauy))[j][i] * I_rho0 * bt.IDatv[j][i] +
          bt.frhatv[0][j][i] * (CAv[0][j][i] + diffv[0][j][i] + PFv[0][j][i]);
      }
      for (k=1;k<=NZ-1;k++) for (i=X1;i<=nx;i++) {
        ubt_force[j][i] += bt.frhatu[k][j][i] * (CAu[k][j][i] + diffu[k][j][i]
                                               + PFu[k][j][i]);
        vbt_force[j][i] += bt.frhatv[k][j][i] * (CAv[k][j][i] + diffv[k][j][i]
                                               + PFv[k][j][i]);
    } }
    btstep(mp,m,dt,ubt_force,vbt_force,u_accel_bt,v_accel_bt);

    for (k=0;k<=NZ-1;k++) for (j=Y1;j<=ny;j++) for (i=X1;i<=nx;i++) {
      v[m][k][j][i] = vmask[j][i] * (v[mp][k][j][i] + dt_pred *
              (PFv[k][j][i] + CAv[k][j][i] + diffv[k][j][i] -
               ((pbce[k][j+1][i]-bt.gtot[j+1][i].S)*bt.etaav[j+1][i] -
                (pbce[k][j][i]-bt.gtot[j][i].N)*bt.etaav[j][i]) * IDYv(j,i) +
               v_accel_bt[j][i]));
    }
    for (k=0;k<=NZ-1;k++) for (j=Y1;j<=ny;j++) for (i=X1;i<=nx;i++) {
      u[m][k][j][i] = umask[j][i] * (u[mp][k][j][i] + dt_pred *
              (PFu[k][j][i] + CAu[k][j][i] + diffu[k][j][i] -
               ((pbce[k][j][i+1]-bt.gtot[j][i+1].W)*bt.etaav[j][i+1] -
                (pbce[k][j][i]-bt.gtot[j][i].E)*bt.etaav[j][i]) * IDXu(j,i) +
               u_accel_bt[j][i]));
    }
    if (calc_bbl) set_viscous_BBL(u[mp],v[mp],h[mp],thermovar);
    calc_bbl = 0;

    vertvisc(u[m],v[m],h[mp],fluxes,bt,dt_pred,u[2],v[2]);

    pass_var(u[m],NZ,1,2,0.0,1,0,TO_NORTH+TO_SOUTH+TO_EAST+TO_WEST,U_PT+VECTOR);
    pass_var(v[m],NZ,1,2,0.0,0,0,TO_NORTH+TO_SOUTH+TO_EAST+TO_WEST,V_PT+VECTOR);
    pass_var(u[2],NZ,1,2,0.0,0,0,TO_NORTH+TO_SOUTH+TO_EAST+TO_WEST,U_PT+VECTOR);
    pass_var(v[2],NZ,1,2,0.0,0,0,TO_NORTH+TO_SOUTH+TO_EAST+TO_WEST,V_PT+VECTOR);

    continuity(u[2],v[2],h[mp],h[m],uh,vh,dt);

    pass_var(h[m],NZ,1,2,EPSILON,1,0,TO_NORTH+TO_SOUTH+TO_EAST+TO_WEST,0);
    pass_var(uh,NZ,1,1,0.0,0,0,TO_SOUTH,U_PT+VECTOR);
    pass_var(vh,NZ,1,1,0.0,0,0,TO_WEST,V_PT+VECTOR);

    for (k=0;k<=NZ-1;k++) for (j=Y1-2;j<=ny+2;j++) for (i=X1-2;i<=nx+2;i++) {
      h[2][k][j][i] = 0.5*(h[0][k][j][i] + h[1][k][j][i]);
    }

/* The correction phase of the time step starts here.                 */
    enable_averaging(dt,rel_time);

    btcalc(u[m],v[m],h[m],m,0);
    if (BEGW != 0.0) {
      for (k=0;k<=NZ-1;k++) for (j=Y1-1;j<=ny+1;j++) for (i=X1-1;i<=nx+1;i++) {
        h[m][k][j][i] = (1.0-BEGW)*h[mp][k][j][i] + (BEGW)*h[m][k][j][i];
      }
      PressureForce(h[m],thermovar,bt,PFu,PFv,pbce,1,0);
    }

    horizontal_viscosity(u[2], v[2], h[2], diffu, diffv);
    CorAdCalc(u[2],v[2],h[2],uh,vh,CAu,CAv);
/* Calculate the momentum forcing terms for the barotropic equations. */
    for (j=Y1;j<=ny;j++) {
      for (i=X1;i<=nx;i++) {
        ubt_force[j][i] = (*(fluxes.taux))[j][i] * I_rho0 * bt.IDatu[j][i] +
          bt.frhatu[0][j][i] * (CAu[0][j][i] + diffu[0][j][i] + PFu[0][j][i]);
        vbt_force[j][i] = (*(fluxes.tauy))[j][i] * I_rho0 * bt.IDatv[j][i] +
          bt.frhatv[0][j][i] * (CAv[0][j][i] + diffv[0][j][i] + PFv[0][j][i]);
      }
      for (k=1;k<=NZ-1;k++) for (i=X1;i<=nx;i++) {
        ubt_force[j][i] += bt.frhatu[k][j][i] * (CAu[k][j][i] + diffu[k][j][i]
                                               + PFu[k][j][i]);
        vbt_force[j][i] += bt.frhatv[k][j][i] * (CAv[k][j][i] + diffv[k][j][i]
                                               + PFv[k][j][i]);
    } }
    btstep(mp,m,dt,ubt_force,vbt_force,u_accel_bt,v_accel_bt);

    for (k=0;k<=NZ-1;k++) for (j=Y1;j<=ny;j++) for (i=X1;i<=nx;i++) {
      double PFu1;
      PFu1 = ((pbce[k][j][i+1]-bt.gtot[j][i+1].W)*bt.etaav[j][i+1] -
              (pbce[k][j][i]-bt.gtot[j][i].E)*bt.etaav[j][i]) * IDXu(j,i);
      u[m][k][j][i] = umask[j][i] * (u[mp][k][j][i] + dt *
            (PFu[k][j][i] + CAu[k][j][i] + diffu[k][j][i] - PFu1 +
             u_accel_bt[j][i]));

      if (diag.PFu_tot != NULL)
        diag.PFu_tot[k][j][i] = PFu[k][j][i] - PFu1 + diag.PFu_bt[j][i];
      if (diag.CAu_tot != NULL)
        diag.CAu_tot[k][j][i] = CAu[k][j][i] + u_accel_bt[j][i]
                             - diag.PFu_bt[j][i];
    }

    for (k=0;k<=NZ-1;k++) for (j=Y1;j<=ny;j++) for (i=X1;i<=nx;i++) {
      double PFv1;
      PFv1 = ((pbce[k][j+1][i]-bt.gtot[j+1][i].S)*bt.etaav[j+1][i] -
              (pbce[k][j][i]-bt.gtot[j][i].N)*bt.etaav[j][i]) * IDYv(j,i);
      v[m][k][j][i] = vmask[j][i] * (v[mp][k][j][i] + dt *
            (PFv[k][j][i] + CAv[k][j][i] + diffv[k][j][i] - PFv1 +
             v_accel_bt[j][i]));

      if (diag.PFv_tot != NULL)
        diag.PFv_tot[k][j][i] = PFv[k][j][i] - PFv1 + diag.PFv_bt[j][i];
      if (diag.CAv_tot != NULL)
        diag.CAv_tot[k][j][i] = CAv[k][j][i] + v_accel_bt[j][i]
                              - diag.PFv_bt[j][i];
    }
    vertvisc(u[m],v[m],h[mp],fluxes,bt,dt,u[2],v[2]);

/*   Here various terms used in to update the momentum equations are  */
/* offered for averaging.                                             */
    {
      static int ref_PFut = -1, ref_PFvt = -1, ref_CAut = -1, ref_CAvt = -1;
      if (ref_PFut < 0) {
        ref_PFut = get_field_ref(&diag.PFu_tot[0]);
        ref_PFvt = get_field_ref(&diag.PFv_tot[0]);
        ref_CAut = get_field_ref(&diag.CAu_tot[0]);
        ref_CAvt = get_field_ref(&diag.CAv_tot[0]);
      }
      if (ref_PFut > 0) average_field(&diag.PFu_tot[0], ref_PFut);
      if (ref_PFvt > 0) average_field(&diag.PFv_tot[0], ref_PFvt);
      if (ref_CAut > 0) average_field(&diag.CAu_tot[0], ref_CAut);
      if (ref_CAvt > 0) average_field(&diag.CAv_tot[0], ref_CAvt);
    }

    pass_var(u[m],NZ,1,2,0.0,1,0,TO_NORTH+TO_SOUTH+TO_EAST+TO_WEST,U_PT+VECTOR);
    pass_var(v[m],NZ,1,2,0.0,0,0,TO_NORTH+TO_SOUTH+TO_EAST+TO_WEST,V_PT+VECTOR);
    pass_var(u[2],NZ,1,2,0.0,0,0,TO_NORTH+TO_SOUTH+TO_EAST+TO_WEST,U_PT+VECTOR);
    pass_var(v[2],NZ,1,2,0.0,0,0,TO_NORTH+TO_SOUTH+TO_EAST+TO_WEST,V_PT+VECTOR);

    continuity(u[2],v[2],h[mp],h[m],uh,vh,dt);

    pass_var(h[m],NZ,1,2,EPSILON,1,0,TO_NORTH+TO_SOUTH+TO_EAST+TO_WEST,0);
    pass_var(uh,NZ,1,1,0.0,0,0,TO_SOUTH,U_PT+VECTOR);
    pass_var(vh,NZ,1,1,0.0,0,0,TO_WEST,V_PT+VECTOR);

    for (k=0;k<=NZ-1;k++) {
      for (j=Y1-2;j<=ny+2;j++) for (i=X1-2;i<=nx+2;i++) {
          h[2][k][j][i] = 0.5*(h[0][k][j][i] + h[1][k][j][i]);
      }
      for (j=Y1-1;j<=ny;j++) for (i=X1-1;i<=nx;i++) {
        uhtr[k][j][i] += uh[k][j][i]; vhtr[k][j][i] += vh[k][j][i];
      }
    }

/*   Here the thickness fluxes are offered for averaging.             */
    {
      static int ref_uh = -1, ref_vh = -1;
      if (ref_uh < 0) {
        ref_uh = get_field_ref(&uh[0]); ref_vh = get_field_ref(&vh[0]);
      }
      if (ref_uh > 0) average_field(&uh[0], ref_uh);
      if (ref_vh > 0) average_field(&vh[0], ref_vh);
    }
    disable_averaging();

#else /* -------------------------------------------------- not SPLIT */

/* Matsuno's third order accurate three step scheme is used to step   */
/* all of the fields except h.  h is stepped separately.              */
    disable_averaging();

    enable_averaging(dt,rel_time);
    horizontal_viscosity(u[0], v[0], h[0], diffu, diffv);
    disable_averaging();

    continuity(u[0],v[0],h[0],h[1],uh,vh,(dt*0.5));
    pass_var(h[1],NZ,1,2,EPSILON,1,0,TO_NORTH+TO_SOUTH+TO_EAST+TO_WEST,0);
    pass_var(uh,NZ,1,1,0.0,0,0,TO_SOUTH,U_PT+VECTOR);
    pass_var(vh,NZ,1,1,0.0,0,0,TO_WEST,V_PT+VECTOR);

    enable_averaging(0.5*dt,rel_time-0.5*dt);
/*   Here the thickness fluxes are offered for averaging.             */
    {
      static int ref_uh = -1, ref_vh = -1;
      if (ref_uh < 0) {
        ref_uh = get_field_ref(&uh[0]); ref_vh = get_field_ref(&vh[0]);
      }
      if (ref_uh > 0) average_field(&uh[0], ref_uh);
      if (ref_vh > 0) average_field(&vh[0], ref_vh);
    }
    disable_averaging();

    for (k=0;k<=NZ-1;k++) {
      for (j=Y1-2;j<=ny+2;j++) for (i=X1-2;i<=nx+2;i++) {
        h[2][k][j][i] = (h[0][k][j][i] + h[1][k][j][i]) * 0.5;
      }
      for (j=Y1;j<=ny;j++) for (i=X1;i<=nx;i++) {
        u[0][k][j][i] += dt * diffu[k][j][i] * umask[j][i];
        v[0][k][j][i] += dt * diffv[k][j][i] * vmask[j][i];
      }
      for (j=Y1-1;j<=ny;j++) for (i=X1-1;i<=nx;i++) {
        uhtr[k][j][i] += 0.5*uh[k][j][i];
        vhtr[k][j][i] += 0.5*vh[k][j][i];
      }
    }
    pass_var(u[0],NZ,1,2,0.0,1,0,TO_NORTH+TO_SOUTH+TO_EAST+TO_WEST,U_PT+VECTOR);
    pass_var(v[0],NZ,1,2,0.0,0,0,TO_NORTH+TO_SOUTH+TO_EAST+TO_WEST,V_PT+VECTOR);

    CorAdCalc(u[0],v[0],h[2],uh,vh,CAu,CAv);
    PressureForce(h[2],thermovar,bt,PFu,PFv,NULL,0,0);
    for (k=0;k<=NZ-1;k++) for (j=Y1;j<=ny;j++) for (i=X1;i<=nx;i++) {
      u[1][k][j][i] = umask[j][i] * (u[0][k][j][i] + dt_pred *
              (PFu[k][j][i] + CAu[k][j][i]));
      v[1][k][j][i] = vmask[j][i] * (v[0][k][j][i] + dt_pred *
              (PFv[k][j][i] + CAv[k][j][i]));
    }
    if (calc_bbl) { set_viscous_BBL(u[0],v[0],h[2],thermovar); calc_bbl = 0;}
    vertvisc(u[1],v[1],h[2],fluxes,bt,dt*0.5,u[2],v[2]);
    pass_var(u[1],NZ,1,2,0.0,1,0,TO_NORTH+TO_SOUTH+TO_EAST+TO_WEST,U_PT+VECTOR);
    pass_var(v[1],NZ,1,2,0.0,0,0,TO_NORTH+TO_SOUTH+TO_EAST+TO_WEST,V_PT+VECTOR);

    continuity(u[1],v[1],h[1],h[2],uh,vh,(0.5*dt));
    pass_var(h[2],NZ,1,2,EPSILON,1,0,TO_NORTH+TO_SOUTH+TO_EAST+TO_WEST,0);
    pass_var(uh,NZ,1,1,0.0,0,0,TO_SOUTH,U_PT+VECTOR);
    pass_var(vh,NZ,1,1,0.0,0,0,TO_WEST,V_PT+VECTOR);

    for (k=0;k<=NZ-1;k++) for (j=Y1-2;j<=ny+2;j++) for (i=X1-2;i<=nx+2;i++) {
      h[2][k][j][i] = (h[1][k][j][i] + h[2][k][j][i]) * 0.5;
    }
    CorAdCalc(u[1],v[1],h[2],uh,vh,CAu,CAv);
    PressureForce(h[2],thermovar,bt,PFu,PFv,NULL,0,0);
    for (k=0;k<=NZ-1;k++) for (j=Y1;j<=ny;j++) for (i=X1;i<=nx;i++) {
      u[2][k][j][i] = umask[j][i] * (u[0][k][j][i] + dt * 0.5 *
              (PFu[k][j][i] + CAu[k][j][i]));
      v[2][k][j][i] = vmask[j][i] * (v[0][k][j][i] + dt * 0.5 *
              (PFv[k][j][i] + CAv[k][j][i]));
    }
    vertvisc(u[2],v[2],h[1],fluxes,bt,dt*0.5,u[2],v[2]);
    pass_var(u[2],NZ,1,2,0.0,1,0,TO_NORTH+TO_SOUTH+TO_EAST+TO_WEST,U_PT+VECTOR);
    pass_var(v[2],NZ,1,2,0.0,0,0,TO_NORTH+TO_SOUTH+TO_EAST+TO_WEST,V_PT+VECTOR);

    continuity(u[2],v[2],h[1],h[0],uh,vh,(dt*0.5));
    pass_var(h[0],NZ,1,2,EPSILON,1,0,TO_NORTH+TO_SOUTH+TO_EAST+TO_WEST,0);
    pass_var(uh,NZ,1,1,0.0,0,0,TO_SOUTH,U_PT+VECTOR);
    pass_var(vh,NZ,1,1,0.0,0,0,TO_WEST,V_PT+VECTOR);

    enable_averaging(0.5*dt,rel_time);
/*   Here the thickness fluxes are offered for averaging.             */
    {
      static int ref_uh = -1, ref_vh = -1;
      if (ref_uh < 0) {
        ref_uh = get_field_ref(&uh[0]); ref_vh = get_field_ref(&vh[0]);
      }
      if (ref_uh > 0) average_field(&uh[0], ref_uh);
      if (ref_vh > 0) average_field(&vh[0], ref_vh);
    }
    disable_averaging();

    for (k=0;k<=NZ-1;k++) {
      for (j=Y1-2;j<=ny+2;j++) for (i=X1-2;i<=nx+2;i++) {
        h[2][k][j][i] = (h[0][k][j][i] + h[1][k][j][i]) * 0.5;
      }
      for (j=Y1-1;j<=ny;j++) for (i=X1-1;i<=nx;i++) {
        uhtr[k][j][i] += 0.5*uh[k][j][i];
        vhtr[k][j][i] += 0.5*vh[k][j][i];
      }
    }

    enable_averaging(dt,rel_time);

    CorAdCalc(u[2],v[2],h[2],uh,vh,CAu,CAv);
    PressureForce(h[2],thermovar,bt,PFu,PFv,NULL,0,0);
    for (k=0;k<=NZ-1;k++) for (j=Y1;j<=ny;j++) for (i=X1;i<=nx;i++) {
      u[0][k][j][i] = umask[j][i] * (u[0][k][j][i] + dt *
              (PFu[k][j][i] + CAu[k][j][i]));
      v[0][k][j][i] = vmask[j][i] * (v[0][k][j][i] + dt *
              (PFv[k][j][i] + CAv[k][j][i]));
    }
    vertvisc(u[0],v[0],h[2],fluxes,bt,dt,u[2],v[2]);
/*   Here various terms used in to update the momentum equations are  */
/* offered for averaging.                                             */
    {
      static int ref_PFu = -1, ref_PFv = -1, ref_CAu = -1, ref_CAv = -1;
      if (ref_PFu < 0) {
        ref_PFu = get_field_ref(&PFu[0]);
        ref_PFv = get_field_ref(&PFv[0]);
        ref_CAu = get_field_ref(&CAu[0]);
        ref_CAv = get_field_ref(&CAv[0]);
      }
      if (ref_PFu > 0) average_field(&PFu[0], ref_PFu);
      if (ref_PFv > 0) average_field(&PFv[0], ref_PFv);
      if (ref_CAu > 0) average_field(&CAu[0], ref_CAu);
      if (ref_CAv > 0) average_field(&CAv[0], ref_CAv);
    }
    disable_averaging();

    pass_var(u[0],NZ,1,2,0.0,1,0,TO_NORTH+TO_SOUTH+TO_EAST+TO_WEST,U_PT+VECTOR);
    pass_var(v[0],NZ,1,2,0.0,0,0,TO_NORTH+TO_SOUTH+TO_EAST+TO_WEST,V_PT+VECTOR);

#endif /* ------------------------------------------------- end SPLIT */

#ifdef THICKNESSDIFFUSE
    thickness_diffuse(h[m],uhtr,vhtr,dt);
    pass_var(h[m],NZ,1,2,EPSILON,1,0,TO_NORTH+TO_SOUTH+TO_EAST+TO_WEST,0);
#endif

    dtnt += dt;
    if ((n%ntstep == 0) || (n==n_max)) {
      enable_averaging(dtnt,rel_time);

      advect_tracer(h[m],uhtr,vhtr,dt);

# ifndef ADIABATIC
      diabatic(u[m],v[m],h[m],thermovar,fluxes,dtnt);

      pass_var(u[m],NZ,1,2,0.0,1,0,TO_NORTH+TO_SOUTH+TO_EAST+TO_WEST,U_PT+VECTOR);
      pass_var(v[m],NZ,1,2,0.0,0,0,TO_NORTH+TO_SOUTH+TO_EAST+TO_WEST,V_PT+VECTOR);
      pass_var(h[m],NZ,1,2,EPSILON,0,0,TO_NORTH+TO_SOUTH+TO_EAST+TO_WEST,0);
#  ifdef BULKMIXEDLAYER
      pass_var(Rml,NML+1,1,1,0.0,0,0,TO_SOUTH+TO_WEST,0);
#  endif
#  ifdef NONLINEAR_EOS
      pass_var(T,NZ,1,1,0.0,0,0,TO_SOUTH+TO_WEST,0);
      pass_var(S,NZ,1,1,0.0,0,0,TO_SOUTH+TO_WEST,0);
#  endif
# endif                                                  /* ADIABATIC */

      calculate_diagnostic_fields(m,dtnt);
      disable_averaging();

      for (k=0;k<=NZ-1;k++) for (i=X1-1;i<=nx;i++) for (j=Y1-1;j<=ny;j++) {
        uhtr[k][j][i] = 0.0; vhtr[k][j][i] = 0.0;
      }
      dtnt = 0.0;
      calc_bbl = 1;
    }

    enable_averaging(dt,rel_time);
/*   Here some of the prognostic variables are offered for averaging. */
    {
      static int ref_u = -1, ref_v = -1, ref_h = -1;
      if (ref_u < 0) {
        ref_u = get_field_ref(&u[0][0]); ref_v = get_field_ref(&v[0][0]);
        ref_h = get_field_ref(&h[0][0]);
      }
      if (ref_u > 0) average_field(&u[m][0], ref_u);
      if (ref_v > 0) average_field(&v[m][0], ref_v);
      if (ref_h > 0) average_field(&h[m][0], ref_h);
    }
    disable_averaging();
  }

  enable_averaging(dt*n_max,rel_time);
  average_forcing(fluxes,dt*n_max);
  disable_averaging();

  calculate_surface_state(u[m], v[m], h[m], state);

  return(m);
}


void HIM_initialize(time_type *day, struct surface *state) {
/* This subroutine takes care of initializing HIM.  The user-         */
/* configurable portions are found in initialize.c, which is called   */
/* from this subroutine.                                              */

/* Arguments: day - Time of the start of the run units of TIMEUNIT.   */
/*  (out)     state - A structure containing fields that describe the */
/*                    initial surface state of the ocean.             */

  int i, j, k;

/*   Initialize the computational domain.                             */
  set_up_parallel();

  thermovar.T = NULL; thermovar.S = NULL; thermovar.Rml = NULL;
  thermovar.frazil = NULL;
  bt.ubtav = NULL; bt.vbtav = NULL; bt.IDatu = NULL; bt.IDatv = NULL;
  bt.gtot = NULL; bt.frhatu = NULL; bt.frhatv = NULL;

/* Set up pointers to all advective tracers.                          */
/*   Zero out the tracer advection velocities.                        */
  for (k=0;k<=NZ-1;k++) for (j=Y1-1;j<=ny;j++) for (i=X1-1;i<=nx;i++) {
    uhtr[k][j][i] = 0.0; vhtr[k][j][i] = 0.0;
  }
#ifdef TEMPERATURE
  register_tracer(&T,NZ,"T"); register_tracer(&S,NZ,"S");
  thermovar.T = &T[0]; thermovar.S = &S[0];
# ifdef FRAZIL
  thermovar.frazil = &frazil[0];
# endif
#endif
#ifdef BULKMIXEDLAYER
  register_tracer(&Rml,NML+1,"Rml");
  thermovar.Rml = &Rml[0];
#endif

/*   Set the fields that are needed for bitwise identical restarting  */
/* the time stepping scheme.                                          */
  set_restart_fields();

/*   Initialize all of the relevant fields.                           */
  initialize(day);


  calculate_surface_state(u[0], v[0], h[0], state);

/*   Update the values at the boundaries of variables that are set    */
/* in initialize.                                                     */
  pass_var(u[0],NZ,1,2,0.0,1,0,TO_NORTH+TO_SOUTH+TO_EAST+TO_WEST,U_PT+VECTOR);
  pass_var(v[0],NZ,1,2,0.0,0,0,TO_NORTH+TO_SOUTH+TO_EAST+TO_WEST,V_PT+VECTOR);
  pass_var(h[0],NZ,1,2,EPSILON,0,0,TO_NORTH+TO_SOUTH+TO_EAST+TO_WEST,0);
#  ifdef BULKMIXEDLAYER
  pass_var(Rml,NML+1,1,1,0.0,0,0,TO_SOUTH+TO_WEST,0);
#  endif
#  ifdef NONLINEAR_EOS
  pass_var(T,NZ,1,1,0.0,0,0,TO_SOUTH+TO_WEST,0);
  pass_var(S,NZ,1,1,0.0,0,0,TO_SOUTH+TO_WEST,0);
#  endif

/*   Calculate several fields which might have been carried in the    */
/* restart file, in cases where a restart file was not used.          */
#ifdef SPLIT                                           /* start SPLIT */
  barotropic_init(u[0], v[0], h[0], &bt);
  if (!query_initialized(diffu[0][0]) || !query_initialized(diffv[0][0]))
     horizontal_viscosity(u[0], v[0], h[0], diffu, diffv);
  if ((!query_initialized(u[2][0][0])) || (!query_initialized(v[2][0][0])))
    for (k=0;k<=NZ-1;k++) for (j=Y1;j<=ny;j++) for (i=X1;i<=nx;i++) {
      u[2][k][j][i] = u[0][k][j][i];
      v[2][k][j][i] = v[0][k][j][i];
    }

/* This call is just here to initialize uh and vh.                    */
  if (!query_initialized(uh[0][0]) || !query_initialized(vh[0][0])) {
    continuity(u[0],v[0],h[0],h[1],uh,vh,(HIM_params.dt));
    pass_var(h[1],NZ,1,2,EPSILON,1,0,TO_NORTH+TO_SOUTH+TO_EAST+TO_WEST,0);
    for (k=0;k<=NZ-1;k++) for (j=Y1;j<=ny;j++) for (i=X1;i<=nx;i++)
      h[2][k][j][i] = 0.5*(h[0][k][j][i] + h[1][k][j][i]);
  } else {
    if (!query_initialized(h[2][0][0]))
      for (k=0;k<=NZ-1;k++) for (j=Y1;j<=ny;j++) for (i=X1;i<=nx;i++)
        h[2][k][j][i] = h[0][k][j][i];
  }

  pass_var(u[2],NZ,1,2,0.0,1,0,TO_NORTH+TO_SOUTH+TO_EAST+TO_WEST,U_PT+VECTOR);
  pass_var(v[2],NZ,1,2,0.0,0,0,TO_NORTH+TO_SOUTH+TO_EAST+TO_WEST,V_PT+VECTOR);
  pass_var(h[2],NZ,1,2,EPSILON,0,0,TO_NORTH+TO_SOUTH+TO_EAST+TO_WEST,0);
  pass_var(uh,NZ,1,1,0.0,1,0,TO_SOUTH+TO_EAST,U_PT+VECTOR);
  pass_var(vh,NZ,1,1,0.0,0,0,TO_WEST+TO_NORTH,V_PT+VECTOR);
#endif                                                   /* end SPLIT */

  calculate_diagnostic_fields(0,0.0);

}

static void set_restart_fields(void) {
/*   Set the fields that are needed for bitwise identical restarting  */
/* the time stepping scheme.  In addition to those specified here     */
/* directly, there may be fields related to the forcing or to the     */
/* barotropic solver that are needed; these are specified in sub-     */
/* routines that are called from this one.                            */
/*   This routine should be altered if there are any changes to the   */
/* time stepping scheme.  The CHECK_RESTART facility may be used to   */
/* confirm that all needed restart fields have been included.         */
  {
    vardesc vd[] = {
       {"u","Zonal velocity",'u','L','s',"meter second-1", 'd'},
       {"v","Meridional velocity",'v','L','s',"meter second-1", 'd'},
       {"h","Layer Thickness",'h','L','s',"meter", 'd'}};
    register_restart_field(u[0][0][0], u[1][0][0], vd[0], 1);
    register_restart_field(v[0][0][0], v[1][0][0], vd[1], 1);
    register_restart_field(h[0][0][0], h[1][0][0], vd[2], 1);
  }
#ifdef SPLIT
  {
    vardesc vd[] = {
       {"u2","Auxiliary Zonal velocity",'u','L','s',"meter second-1", 'd'},
       {"v2","Auxiliary Meridional velocity",'v','L','s',"meter second-1", 'd'},
       {"h2","Auxiliary Layer Thickness",'h','L','s',"meter", 'd'},
       {"uh","Zonal thickness flux",'u','L','s',"meter3 second-1", 'd'},
       {"vh","Meridional thickness flux",'v','L','s',"meter3 second-1", 'd'},
       {"diffu","Zonal horizontal viscous acceleration",'u','L','s',
           "meter second-2", 'd'},
       {"diffv","Meridional horizontal viscous acceleration",'v','L','s',
           "meter second-2", 'd'}};
    register_restart_field(u[2][0][0], NULL, vd[0], 0);
    register_restart_field(v[2][0][0], NULL, vd[1], 0);
    register_restart_field(h[2][0][0], NULL, vd[2], 0);
    register_restart_field(uh[0][0], NULL, vd[3], 0);
    register_restart_field(vh[0][0], NULL, vd[4], 0);
    register_restart_field(diffu[0][0], NULL, vd[5], 0);
    register_restart_field(diffv[0][0], NULL, vd[6], 0);
    register_barotropic_restarts();
  }
#endif
#ifdef BULKMIXEDLAYER
  { vardesc vd = {"Rml","Mixed Layer Potential Density",'h','2','s',
                  "kg meter-3", 'd'};
    register_restart_field(Rml[0][0], NULL, vd, 1); }
#endif
#ifdef TEMPERATURE
  { vardesc vdT = {"Temp","Potential Temperature",'h','L','s',"degC", 'd'};
    vardesc vdS = {"Salt","Salinity",'h','L','s',"PSU", 'd'};
    register_restart_field(T[0][0], NULL, vdT, 1);
    register_restart_field(S[0][0], NULL, vdS, 1); }
#endif

  register_forcing_restarts();
}

void calculate_surface_state(double u[NZ][NYMEM][NXMEM],
         double v[NZ][NYMEM][NXMEM], double h[NZ][NYMEM][NXMEM],
         struct surface *state) {
/*   This subroutine sets the surface (return) properties of the ocean*/
/* model by setting the appropriate pointers in state.  Unused fields */
/* are set to NULL.                                                   */

/* Arguments: u - Zonal velocity, in m s-1.                           */
/*  (in)      v - Meridional velocity, in m s-1.                      */
/*  (in)      h - Layer thickness, in m.                              */
/*  (out)     state - A structure containing fields that describe the */
/*                    surface state of the ocean.                     */

#ifndef BULKMIXEDLAYER
  double depth[NXMEM];        /* The distance from the surface, in m. */
  double depth_ml;            /* The depth over which to average to   */
                              /* determine mixed layer properties, m. */
                              /* their climatological values.         */
  static double Tml[NYMEM][NXMEM];/* The temperature averaged over the*/
                              /* fluid that is within a mixed layer   */
                              /* thickness of the surface, in C.      */
  static double Sml[NYMEM][NXMEM];/* The salinity averaged over the   */
                              /* fluid that is within a mixed layer   */
                              /* thickness of the surface, in PSU.    */
  static double Rhoml[NYMEM][NXMEM];/* The potential density averaged */
                            /* over the fluid that is within a mixed  */
                            /* layer thickness of the surface, kg m-3.*/
  int i, j, k;

  depth_ml = HIM_params.Hmix;

/*   Determine the mean properties of the uppermost depth_ml fluid.   */
  for (j=Y1;j<=ny;j++) {
    for (i=X1;i<=nx;i++) {
      depth[i] = 0.0;
      Tml[j][i] = 0.0;
      Sml[j][i] = 0.0;
      Rhoml[j][i] = 0.0;
    }
    for (k=0;k<=NZ-1;k++) {
      for (i=X1;i<=nx;i++) {
        double dh;  /* The thickness of a layer that is within the    */
                    /* mixed layer, in m.                             */
        if (depth[i]+h[k][j][i] < depth_ml) dh = h[k][j][i];
        else if (depth[i] < depth_ml) dh = depth_ml-depth[i];
        else dh = 0.0;
        if (thermovar.T != NULL) {
          Tml[j][i] += dh * thermovar.T[k][j][i];
          Sml[j][i] += dh * thermovar.S[k][j][i];
        }
        else Rhoml[j][i] += dh * Rlay[k];
        depth[i] += dh;
      }
    }
/* Calculate the average properties of the mixed layer depth.         */
    for (i=X1;i<=nx;i++) {
      if (thermovar.T != NULL) {
        Tml[j][i] /= depth[i];
        Sml[j][i] /= depth[i];
      }
      else Rhoml[j][i] /= depth[i];
    }
  }

/* Any unused elements of state are set to NULL.                      */
  if (thermovar.T != NULL) {
    state->SST = &Tml;
    state->SSS = &Sml;
    state->Rml = NULL;
  } else {
    state->SST = NULL;
    state->SSS = NULL;
    state->Rml = &Rhoml;
  }

#else                                       /* end NOT BULKMIXEDLAYER */
/* Any unused elements of state are set to NULL.                      */
  if (thermovar.T != NULL) {
    state->SST = &thermovar.T[0];
    state->SSS = &thermovar.S[0];
    state->Rml = NULL;
  } else {
    state->SST = NULL;
    state->SSS = NULL;
    state->Rml = &thermovar.Rml[0];
  }
#endif                                          /* end BULKMIXEDLAYER */

  state->u = &u[0];
  state->v = &v[0];
  state->frazil = (double (*)[NYMEM][NXMEM]) thermovar.frazil;
}
