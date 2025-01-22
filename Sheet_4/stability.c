#include "grid/multigrid.h"
#include "embed.h"
#include "navier-stokes/centered.h"
#include "view.h"
#include "vof.h"
#include "two-phase.h"
#include "tension.h"
#include "fractions.h"

#include <omp.h>

#define EPS 1e-14
#define maxlevel 5
#define Pi 3.14159265359


double Reynolds = 100.;
// !RC Specifying an end time is helpful
double t_end = 5.01;
double period = 2;

face vector muv[];

// !RC Runtime statistics output stream
FILE * fp_stats;
FILE * fp_stats_2;

// !RC Boundary condition setup worth doing pen and paper first on a schematic
u.t[top] = dirichlet(1.);
u.n[top] = dirichlet(0.);
u.t[embed] = dirichlet(0.);
u.n[embed] = dirichlet(0.);

int main(int argc, char **argv)
{
  
  size (8. [0]);
  origin (-L0/2.0, 0.5-L0);
  //N = atoi(argv[1]); // sets grid points per dimension
  N = 128;
  fprintf(stderr,"%d",N);
  mu = muv; // what does this line do??
  periodic(right);
  
  // !RC Pointer of the file to save stats
  {
    char name[200];
    sprintf(name, "logstats.dat");
    fp_stats = fopen(name, "w");
  } 

   // information at intervals
  {
    char name_2[200];
    sprintf(name_2, "VerticleVelocity.dat");
    fp_stats_2 = fopen(name_2, "w");
  } 
  dt = 0.001;
  run(); 
  
  fclose(fp_stats);
  fclose(fp_stats_2);
}

event properties (i++)
{
  foreach_face()
    // !RC The viscosity can be thought of as an inverse Reynolds number
    muv.x[] = fm.x[]/Reynolds; 
}



scalar un[];

scalar u_x_previous[];

event init (t = 0.) {
  // Use embedded boundary to construct channel geometry
  // !RC A lesser known "feature" of the code is that embedded geometries shouldn't be aligned with natural gridcells, 
  // hence the 0.501 instead of 0.5 here (will explain more at our next meeting)
  solid (cs, fs, EPS+0.5 + y);

  //defining the interface
  fraction(f,0.1*sin(period*Pi*2*x/8)+0.1-y);

  //defining densities and viscosity
  //rho1 = 1;
  mu1 = 1;
  //rho2 = 2;
  mu2 = 2;

  foreach(){
    	u.x[] = cs[] ? 0. : 0.;
    	un[]=u.x[];
  }
  
}


 // Check for convergence
event logfile (i++) {
  // change(s,sn) computes the max absolute diff between s and sn, then sets sn to s
  double du = change (u.x,un);
  fprintf (stderr, "%d %g %d %d %g\n", i, t, mgp.i, mgu.i,du);
  // checks if a steady state has been reached
  
  if (i>0 && du<1e-2){
    return 1;
  }
  foreach(){
    	un[]=u.x[];
  }
}


// !RC Adaptivity not needed yet
/*
event adapt (i++) {
  adapt_wavelet ({cs,u,f}, (double[]){1e-2,3e-2,3e-2,3e-2}, maxlevel, 4);
}
*/

// !RC Exercise: using file input/output, save the same information at regular intervals, not just at the end
event data_out_final (t = end){
  
    int j;
    double u_x = 0;
    double u_x_2 = 0;
    double y_step = -0.5;

    for (j=0;  j<100; j++) {
      u_x = interpolate(u.x, 0.0, y_step);
      u_x_2 = interpolate(u.x,y_step*8.0,0.0);

      // And write them to file
      printf ( "%g %g \n", y_step, u_x);
      fprintf(fp_stats_2,"%g %g \n", y_step*8.0, u_x_2);
      
      y_step += (1./100.);
      
    }
    // final value slighly bellow
    u_x = interpolate(u.x, 0.0, 0.5-0.001);
    printf("%g %g \n", 0.5, u_x); 


    //show which fluid is witch

    view (fov = 27.0, tx = 0.0, ty = 0.4, width = 2000, height = 2000);
    clear();
    squares ("f", map = cool_warm, min = 0.0, max = 1.0);
    save ("fluids.png");

    
    
}

// saves horizontal speed at x=0 for every 0.5 seconds

event data_out (t += 0.1){
  
    int j;
    double u_x = 0;
    double y_step = -0.5;
    char str[80];

    sprintf(str, "time=%.3f.dat", t);
    FILE * fp = fopen (str, "w");

    if (fp == NULL) {
      fprintf(stderr,"Error opening file!\n");
      exit(1);
    }

    for (j=0;  j<100; j++) {
      u_x = interpolate(u.x, 0.0, y_step);
      // And write them to file
      fprintf (fp, "%g %g \n", y_step, u_x);
      
      y_step += (1./100.);
      
    } 
    fclose(fp);
}

// !RC Regular output stream of runtime statistics
event logstats (t += 0.1) {

    timing s = timer_timing (perf.gt, i, perf.tnc, NULL);
 
    // i, timestep, no of cells, real time elapsed, cpu time
    fprintf(fp_stats, "i: %i t: %g dt: %g #Cells: %ld Wall clock time (s): %g CPU time (s): %g \n", i, t, dt, grid->n, perf.t, s.cpu);
    fflush(fp_stats);
}





// !RC Output animations
event movies (t += 0.05)
{
  char timestring[100];

  view (fov = 27.0, tx = 0.0, ty = 0.4, width = 2000, height = 2000);
  clear();
  squares ("u.x", map = cool_warm, min = 0.0, max = 1.0);
  squares ("cs", map = gray);
  sprintf(timestring, "t=%2.02f",t);
  draw_string(timestring, pos=1, lc= { 0, 0, 0 }, lw=2);
  save ("HorizontalVelocity.mp4");

  view (fov = 27.0, tx = 0.0, ty = 0.4, width = 2000, height = 2000);
  clear();
  squares ("f", map = cool_warm, min = 0.0, max = 1.0);
  squares ("fs", map = gray);
  sprintf(timestring, "t=%2.02f",t);
  draw_string(timestring, pos=1, lc= { 0, 0, 0 }, lw=2);
  save ("Interface.mp4");

}

event failure(t = t_end){
  fprintf (stderr, "Failed to converge to steady state\n");
}