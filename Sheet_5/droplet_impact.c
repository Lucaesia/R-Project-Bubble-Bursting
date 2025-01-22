/* droplet_impact.c
    A Basilisk script to model the impact of a droplet of water impacting onto a
    rigid surface. 

    The problem involves a droplet of diameter 1mm travelling downwards at speed
    1m/s (falling under gravity), surrounded air in an ambient state impacting 
    onto a rigid surface. The contact angle of the droplet on the surface is 
    90°. 

    The problem is non-dimensionalised with regard to the droplet phase 
    parameters i.e. the water phase will have unity density, viscosity etc. 
    Therefore, the parameters used in this simulation are:

    Water density = 1000 kg/m^3, Air density = 1.23 kg/m^3 
        => Density ratio RHO_R = 1.23e-3
    Water viscosity = 8.9e-4 Pa s, Air viscosity = 1.81e-5 Pa s
        => Viscosity ratio MU_R = 2.03e-2
    Surface tension sigma = 72.9e-3 N/m

    Hence the dimensionless numbers which categorise the flow are:
    Reynolds number RE = rho * V0 * D0 / mu1 = 1120
    Weber number WE = rho * V0^2 * D0 / sigma = 13.7
    Froude number FR = V0 / sqrt(g D0) = 1 / sqrt(9.81 * 1e-3) = 10.1
    

    With these parameters, we will have a droplet with dimensionless diameter 1,
    dimensionless speed 1, and an acceleration due to gravity 1 / FR^2. 

    We will assume the flow is axisymmetric. In Basilisk, using axi.h, the x 
    coordinate is the height above the surface z and the y coordinate is the 
    radial distance from the z axis r. Therefore, our droplet is falling in the
    -x direction and spreading in the positive y direction. 
*/

#include "axi.h"
#include "navier-stokes/centered.h"
#include "two-phase.h"
#include "tension.h" 
#include "vof.h"
#include "contact.h"
#include "view.h"
#include "tag.h"
#include <omp.h>

#define MINLEVEL 4
#define MAXLEVEL 11

/* Physical constants */
double RHO_R = 1.23e-3; // Density ratio
double MU_R = 2.03e-2; // Viscosity ratio
double RE = 112.;
double WE = 1370;
double FR = 10.1;
double CONTACT_ANGLE = 90.;
double CENTRE = 0.7;

/* Computational constants */
double BOX_WIDTH = 5.; // Width of the computational box
int output_no = 1; // Number of output files for the interface
int gfs_output_no = 1;

// Parameters for the refined box around the origin
double REFINED_BOX_WIDTH = 0.1; // Width of refined box around origin
double REFINED_BOX_HEIGHT = 0.05; // Height of refined box around origin

// Also refine near the surface where the droplet is expected to spread
double REFINED_SURFACE_HEIGHT = 0.01; 

// Variables to measure the CPU time of the simulation
double  start_time, end_time, total_time;

/* Names of data files */
char general_interface_filename[] = "interface_positions";
char interface_time_filename[] = "interface_times.txt";

/* Define fields */
vector h[]; // Height field of droplet on surface


/* Boundary conditions */

// Conditions for entry from above
u.n[right] = neumann(0.); // Allows outflow through boundary
p[right] = dirichlet(0.); // 0 pressure far from surface

// Conditions on surface
u.n[left] = dirichlet(0.); // No flow through surface
h.t[left] = contact_angle (CONTACT_ANGLE * pi / 180.); // Contact angle
u.t[left] = dirichlet(0.); // No slip at surface

// Conditions far from the droplet in the radial direction
u.n[top] = neumann(0.); // Allows outflow through boundary
p[top] = dirichlet(0.); // 0 pressure far from surface



int main() {
    init_grid(1 << MAXLEVEL); // Create grid according the the level

    size(BOX_WIDTH); // Size of the domain

    /* Set physical constants */
    rho1 = 1.; // Density of water phase
    rho2 = RHO_R; // Density of air phase
    mu1 = 1 / RE; // Viscosity of water phase
    mu2 = mu1 * MU_R; // Viscosity of air phase
    f.sigma = 1./WE; // Surface tension at interface

    run(); // Runs the simulation
}

event init(t = 0) {
/* Initialises the flow as a spherical droplet falling downwards */
    start_time = 0;//omp_get_wtime();; // Time at start of simulation

    // The drop is initially centred at z = 0.6
    fraction(f, - sq(x - CENTRE) - sq(y) + sq(0.5));

    // Associate the height function with the VOF tracer
    f.height = h;

    // The initial velocity is -1 in the z direction
    foreach() {
        u.x[] = -f[];
    }

    // Refine the level around the origin where an entrapped bubble is expected
    refine(x < REFINED_BOX_HEIGHT && y < 0.5 * REFINED_BOX_WIDTH \
        && level < MAXLEVEL);
    
    // Creates the file to save the time at which the interface is outputted
    FILE *interface_time_file = fopen(interface_time_filename, "w");
    fclose(interface_time_file);

}

event logfile (i += 1) {
/* Prints useful quantities to the log file */
    // Prints the timestep number, size and time to the log file
    fprintf(stderr, "i = %d, dt = %g, t = %g \n", i, dt, t );
}

event acceleration (i += 1) {
/* Adds acceleration due to gravity at each time step */
    face vector av = a; // Acceleration at each face
    foreach_face(x) av.x[] -= 1./sq(FR); // Adds acceleration due to gravity
}

event small_droplet_removal (i += 1) {
/* Removes any small droplets that have formed, that are smaller than a specific    
    size */
    remove_droplets(f, 3); // Removes droplets of diameter 3 cells or less
}

event adapt (i++) {
/* Adapts the quadtree grid depending on the error for velocity, vorticity and 
    location of the interface */
    scalar omega[]; // Scalar field for vorticity
	vorticity(u, omega);

    // Adapt the mesh for rapid changes in the free surface, velocity and 
    // vorticity.
    adapt_wavelet ({f, u.x, u.y, omega}, (double[]){1e-4, 5e-2, 5e-2, 1e0},
        maxlevel = MAXLEVEL, minlevel = MINLEVEL);
    
    // Refine the level around the origin where an entrapped bubble is expected
    refine(x < REFINED_BOX_HEIGHT && y < 0.5 * REFINED_BOX_WIDTH \
        && level < MAXLEVEL);

    // Refine the surface where we will have a liquid film
    // refine(x < REFINED_SURFACE_HEIGHT && level < MAXLEVEL);
}

event viewing (t += 1e-3) {
/* Event for creating movies of the evolution of the relevant quantities */

    // Creates a string with the time to put on the plots
    char time_str[80];
    sprintf(time_str, "t = %g\n", t);

    // Sets up box for viewing
    view (width = 1000, height = 1000, fov = 20, ty = -0.5, quat = {0, 0, -0.707, 0.707});

    // Velocity in the vertical direction
    clear();
    draw_vof ("f", lw = 2); // Draws the interface
    squares ("u.x", linear = true); // Linear interpolated vertical velocity
    // box (notics = true);
    mirror ({0,1}) {
        draw_vof ("f", lw = 2);
        squares (color = "u.x", linear = true);
        // box (notics = true);
    }
    draw_string(time_str, pos=1, lc= { 0, 0, 0 }, lw=2);
    save ("vertical_velocity.mp4");

    // Velocity in the horizontal direction
    clear();
    draw_vof ("f", lw = 2); // Draws the interface
    squares ("u.y", linear = true); // Linear interpolated horizontal velocity
    // box (notics = true);
    mirror ({0,1}) {
        draw_vof ("f", lw = 2);
        squares (color = "u.y", linear = true);
        // box (notics = true);
    }
    draw_string(time_str, pos=1, lc= { 0, 0, 0 }, lw=2);
    save ("horizontal_velocity.mp4");

    // Pressure
    clear();
    draw_vof ("f", lw = 2); // Draws the interface
    squares ("p", linear = true); // Linear interpolated pressure
    // box (notics = true);
    mirror ({0,1}) {
        draw_vof ("f", lw = 2);
        squares (color = "p", linear = true);
        // box (notics = true);
    }
    draw_string(time_str, pos=1, lc= { 0, 0, 0 }, lw=2);
    save ("pressure.mp4");

    // Level
    scalar l[];
    foreach() l[] = level;
    clear();
    draw_vof("f", lw = 2);
    squares("l");
    draw_string(time_str, pos=1, lc= { 0, 0, 0 }, lw=2);
    save ("levels.mp4");

    // Colour function
    clear();
    draw_vof("f", lw = 2);
    squares("f", linear = true);
    mirror ({0,1}) {
        draw_vof ("f", lw = 2);
        squares (color = "f", linear = true);
        // box (notics = true);
    }
    draw_string(time_str, pos=1, lc= { 0, 0, 0 }, lw=2);
    save ("vof_tracer.mp4");

}

event interface_calculation (t += 5e-3) {
/* Calculates the relevant quantities relating to the interface position */

    // Creates text file to save output to
    char interface_filename[80];
    sprintf(interface_filename, "%s_%d.txt", \
        general_interface_filename, output_no);
    FILE *interface_file = fopen(interface_filename, "w");

    /* Appends the timestep number and the time for this specific output_no, so 
    that in post-processing we know what these are for each interface file */
    FILE *interface_time_file = fopen(interface_time_filename, "a");
    fprintf(interface_time_file, "%d %d %g\n", output_no, i, t); 
    fclose(interface_time_file);

    // Saves the locations of the interface to the file
    output_facets(f, interface_file);

    output_no++; // Increments the number of outputs there have been

}

event gfsview (t+=1e-2) {
/* Outputs the status of the simulation as a gfs file in order to view later */
    char gfs_filename[80];
    sprintf(gfs_filename, "gfs_output_%d.gfs", gfs_output_no);
    output_gfs(file = gfs_filename);
    gfs_output_no++;
}

event end(t = 50.) {
/* Ends the simulation at a specific time */
    end_time = 1;//omp_get_wtime();; // Time at end of the simulation
    total_time = end_time - start_time; 
    fprintf(stderr, "Total CPU time = %g\n", total_time);
    return 1;
}