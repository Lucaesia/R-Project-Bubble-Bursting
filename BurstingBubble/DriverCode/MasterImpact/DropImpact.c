// Author: Lyllian Chanerley
// Based on code by: Radu Cimpeanu
// Date: 26/09/2022

#include "axi.h"                     // axisymmetric geometry
#include "navier-stokes/centered.h"  // solve NS equations
#define FILTERED                     // Smear density and viscosity jumps
#define mu(f)  (1./(clamp(f,0,1)*(1./mu1 - 1./mu2) + 1./mu2))
#include "two-phase.h"
#include "tension.h"                 // include surface tension between phases
#include "vof.h"                     // solve using VoF method
#include "fractions.h"               // initially defined fraction of fluid phases
#include "view.h"                    // need to make the animations
#include "draw.h"                    // visualisation helper
#include "tag.h"                     // helps track droplet properties
#include <omp.h>
#include <math.h>

               // adds an interpolator

// Dimensional quantities (to be passed as arguments in main below)
double rhoLiquid; // liquid phase density (kg/m^3)
double rhoGas;    // gas phase density (kg/m^3)

double muLiquid;  // liquid dynamic viscosity (kg/ms)
double muGas;     // gas dynamic viscosity(kg/ms)

double sig;       // surface tension (N/m)

double g_accel;   // gravitational acceleration (m/s^2)

double dRadius;   // drop radius (m)

double v_init;    // initial drop velocity (m/s)

// dimensionless specifications (key groupings defined in main below)
#define rho_ratio   (rhoGas/rhoLiquid) // density ratio
#define mu_ratio    (muGas/muLiquid)   // viscosity ratio 

#define poolHeight 10.0                // Pool height (in radii)
#define domainSize 20.0                // Computational box size (in radii)

face vector av[];

FILE * fp_stats;
FILE * fp_vol;
FILE * fp_droplets;
FILE * fp_jet_vel;
FILE * fp_jet_eject;
FILE * menisc_file;
FILE * over_file;
FILE * under_file;
FILE * RL_file;
FILE * fp_interface;


double ND_Weber;
double ND_Reynolds;
double ND_Froude;
double ND_Bond;
double ND_Ohnesorge;
double ND_Laplace;

double filmHeight;

double max_height = -0.901;
int num_jets = 0;

int minLevel = 6;
int maxLevel; // = 11;

double Theta;

double tEnd;


int IC;

// Bottom of Pool = LEFT of domain
// No slip, no permeability
u.n[left] = dirichlet(0.);
u.t[left] = dirichlet(0.);
p[left] = neumann(0.);
pf[left] = neumann(0.);

// Side of Pool = TOP of domain
// Currently no slip and no permeability
// May change to slip condition
u.n[top] = dirichlet(0.); // Impermeability
u.t[top] = neumann(0.); // Slip

// Above the Pool = RIGHT of domain
// Outflow
u.n[right] = neumann(0.);
p[right] = dirichlet(0.);
pf[right] = dirichlet(0.);

// Default for bottom is symmetry

int main(int argc, char * argv[]) {

  rhoLiquid = atof(argv[1]);  // prescribed liquid density 
  rhoGas = atof(argv[2]);     // prescribed gas density
  muLiquid = atof(argv[3]);   // prescribed liquid (dynamic) viscosity
  muGas = atof(argv[4]);      // prescribed gas (dynamic) viscosity
  sig = atof(argv[5]);        // prescribed surface tension coefficient
  g_accel = atof(argv[6]);    // prescribed gravitational acceleration
  dRadius = atof(argv[7]);    // prescribed drop radius
  v_init = atof(argv[8]);     // prescribed initial velocity
  tEnd = atof(argv[9]);       // prescribed simulation end time
  maxLevel = atof(argv[10]);  // prescribed maximum resolution level
  Theta = atof(argv[11]);
  IC = atof(argv[12]);
  
  
  ND_Weber = (rhoLiquid*pow(v_init,2.0)*dRadius)/sig;
  ND_Reynolds = (rhoLiquid*v_init*dRadius)/muLiquid;
  ND_Froude = v_init/pow(dRadius*g_accel,0.5);
  ND_Laplace = ND_Reynolds*ND_Reynolds/ND_Weber;
  ND_Bond = rhoLiquid*g_accel*pow(dRadius,2.0)/sig;
  ND_Ohnesorge = muLiquid/pow(rhoLiquid*sig*dRadius,0.5);
  printf("%d\n",1 << 10)  ;
  init_grid(1 << 9);
  
  size(domainSize);                     
  origin(-0.5*domainSize, 0.0);

  // Create output folders
  mkdir("Slices", 0700);
  mkdir("Animations", 0700);
  mkdir("Interfaces", 0700);

  // Print dimensionless numbers for verification
  //fprintf(stdout, "Reynolds number = %0.6f \n", ND_Reynolds); fflush(stdout);
  //fprintf(stdout, "Weber number = %0.6f \n", ND_Weber); fflush(stdout);
  //fprintf(stdout, "Froude number = %0.6f \n", ND_Froude); fflush(stdout);
  fprintf(stdout, "Bond number = %0.6f \n", ND_Bond); fflush(stdout);
  //fprintf(stdout, "Ohnesorge number = %0.6f \n", ND_Ohnesorge); fflush(stdout);
  fprintf(stdout, "Laplace number = %0.6f \n", ND_Laplace); fflush(stdout);

  rho1 = 1.;
  rho2 = rho_ratio;//gas/liquid
  
  mu1 = 1./pow(ND_Laplace*ND_Bond,0.5);
  mu2 = mu_ratio*mu1;
  
  f.sigma = 1./ND_Bond;
  //f2.sigma = 1./ND_Weber;

  a = av;

  // Pointer of the file to save stats
  {
    char name[200];
    sprintf(name, "logstats.dat");
    fp_stats = fopen(name, "w");
  }

  // Pointer of the file to save interface_stats
  {
    char name[200];
    sprintf(name, "interfacestats.dat");
    fp_interface = fopen(name, "w");
  }

  // Pointer of the file to droplet count info
  {
    char name[200];
    sprintf(name, "logdroplets.dat");
    fp_droplets = fopen(name, "w");
  }

  // Pointer of the file to jet velocity
  {
    char name[200];
    sprintf(name, "jet_vel.dat");
    fp_jet_vel = fopen(name, "w");
  }

  // Pointer to jet ejection data
  {
    char name[200];
    sprintf(name, "jet_eject.dat");
    fp_jet_eject = fopen(name, "w");
  }
  /*--------------------------------*/
  /* Pointers to Initial conditions */
  /*--------------------------------*/
  {
    char name[200];
    sprintf(name, "../../../InitialConditions/theta%.1f/menisc.dat",Theta);
    menisc_file = fopen(name, "r");
    sprintf(name, "../../../InitialConditions/theta%.1f/over.dat",Theta);
    over_file = fopen(name, "r");
    sprintf(name, "../../../InitialConditions/theta%.1f/under.dat",Theta);
    under_file = fopen(name, "r");
    sprintf(name, "../../../InitialConditions/theta%.1f/RL.dat",Theta);
    RL_file = fopen(name, "r");
  }
  

  DT = 1e-2;
  NITERMIN = 1; // default 1
  NITERMAX = 200; // default 100
  TOLERANCE = 1e-6; // default 1e-3
  
  run();

  fclose(fp_stats);
  fclose(fp_droplets);
  fclose(fp_jet_vel);
  fclose(fp_jet_eject);
  fclose(fp_interface);
}

int count_lines(FILE* file)
{
    char buf[65536];
    int counter = 0;
    for(;;)
    {
        size_t res = fread(buf, 1, 65536, file);
        if (ferror(file))
            return -1;

        int i;
        for(i = 0; i < res; i++)
            if (buf[i] == '\n')
                counter++;

        if (feof(file))
            break;
    }

    return counter;
}

int inside_or_out(double x, double y, double x_bar, double x_R, double menisc_max, double* under_x, double * over_x,
                  double* menisc_x, double* under_y, double * over_y, double* menisc_y,
                int under_N,  int over_N, int menisc_N) {
    double val;
    double valp;  // no idea what this does but the spline function want it
    double val2;
    double val3;
    // CONSIDER CASE WHERE MENISC RUNS OUT, done
    //printf("%g %g %g\n", x_bar, x_R, menisc_max);
    if (x>menisc_max){
      if (y>0){
        return 0;
      }
      return 1;
    }

    if ((x > x_bar) && (x > x_R)){
      
      spline_linear_val( menisc_N, menisc_x, menisc_y, x, &val, &valp );
      if (y>val){
        return 0;
      }
      return 1;
    }

    if ((x > x_bar) && (x <= x_R)){
      
      spline_linear_val( menisc_N, menisc_x, menisc_y, x, &val, &valp );
      spline_linear_val( over_N, over_x, over_y, x, &val2, &valp );
      spline_linear_val( under_N, under_x, under_y, x, &val3, &valp );
      if ((y<val)&&(y>val2)){
        return 1;
      }
      if (y<val3){
        return 1;
      }
      return 0;
    }
    // FOR INITIAL AIR BUBBLE
  
    spline_linear_val( under_N, under_x, under_y, x, &val, &valp );
    if (y<val){
      return 1;
    }
    return 0;


}

int inside_or_out_grid_refinement(double x, double y, double x_bar, double x_R, double menisc_max, double* under_x, double * over_x,
  double* menisc_x, double* under_y, double * over_y, double* menisc_y,
  int under_N,  int over_N, int menisc_N) {
  double val;
  double valp;  // no idea what this does but the spline function want it
  double val2;
  double val3;
  // CONSIDER CASE WHERE MENISC RUNS OUT, done
  //printf("%g %g %g\n", x_bar, x_R, menisc_max);
  if (x>menisc_max){
    if (fabs(y)<0.05){
      return 1;
    }
    return 0;
  }

  if ((x > x_bar) && (x > x_R)){

  spline_linear_val( menisc_N, menisc_x, menisc_y, x, &val, &valp );
    if (fabs(y-val)<0.5){
      return 1;
    }
    return 0;
  }

  if ((x > x_bar) && (x <= x_R)){

    spline_linear_val( menisc_N, menisc_x, menisc_y, x, &val, &valp );
    spline_linear_val( over_N, over_x, over_y, x, &val2, &valp );
    spline_linear_val( under_N, under_x, under_y, x, &val3, &valp );
    if ((fabs(y-val)<0.5)||(fabs(y-val2)<0.5)||(fabs(y-val3)<0.5)){
      return 1;
    }
    
    return 0;
  }
  // FOR INITIAL AIR BUBBLE

  spline_linear_val( under_N, under_x, under_y, x, &val, &valp );
  if (fabs(y-val)<0.5){
    return 1;
  }
  return 0;


}



event acceleration (i++) {
  foreach_face(x)  
    av.x[] -= 1;//1./pow(ND_Froude,2.0);
  foreach_face(y)  
    av.y[] += 0.0;
}

event ejection (i++) {
  int state = 0;
  int saved_state = 0;
  int changes = 0;
  double gap_num = 1000;
  double gap = domainSize/gap_num;
  for (double j=-domainSize/2 + 0.000001;j<domainSize/2; j+= gap){
    if (interpolate(f, j,0)<0.5){
      state = 1;
    }
    else{
      state = 0;
    }
    if (state != saved_state){
      changes += 1;
      saved_state = state;
    }
  }
  if (1 == saved_state){
    changes += 1;
    saved_state = state;
  }
  changes = changes/2;
  if (num_jets < changes){
    fprintf(fp_jet_eject,"%d Droplets seperated at time: %g\n", changes-num_jets,t);
    num_jets = changes;
  }
  else if( num_jets > changes){
    fprintf(fp_jet_eject,"%d Droplets disappeared at time: %g\n", num_jets-changes,t);
    num_jets = changes;
  }

}

scalar omega[], viewingfield[], mylevel[], velnorm[];

event init (t = 0.0) {
  /*---------------------------------------------------------------*/
  /* Reading initial condition files and setting initial condition */
  /*---------------------------------------------------------------*/
  double x_bar;
  double x_R;
  double menisc_max;
  int under_N;
  int over_N;
  int menisc_N;
  

  fscanf(RL_file,"%lf\n%lf\n%lf\n%ld\n%ld\n%ld\n", &x_R, &menisc_max, &x_bar, &under_N, &over_N, &menisc_N);//x_bar, x_R, menisc_max
  //printf("hope in the workplace : %d,%d,%d\n ",under_N, over_N, menisc_N);

  double* under_x = malloc(sizeof(double)*(under_N));
  double * over_x = malloc(sizeof(double)*(over_N));
  double* menisc_x = malloc(sizeof(double)*(menisc_N));
  double* under_y = malloc(sizeof(double)*(under_N));
  double * over_y = malloc(sizeof(double)*(over_N));
  double* menisc_y = malloc(sizeof(double)*(menisc_N));
  
  for (int i = 0; i < under_N; i++) {
    double x,y;
    fscanf(under_file,"%lf %lf\n", &x, &y);
    //printf("%g, %g \n",x,y);
    under_x[i] = x;
    under_y[i] = y;
    
  }
  for (int i = 0; i < over_N; i++) {
    double x,y;
    fscanf(over_file,"%lf %lf\n", &x, &y);
    over_x[i] = x;
    over_y[i] = y;
  }
  for (int i = 0; i < menisc_N; i++) {
    double x,y;
    fscanf(menisc_file,"%lf %lf\n", &x, &y);
    menisc_x[i] = x;
    menisc_y[i] = y;
    
  }
  
  
  fclose(menisc_file);
  fclose(over_file);
  fclose(under_file);
  fclose(RL_file);

  filmHeight = -domainSize/2. + poolHeight;
  printf("IC = %d\n",IC);
  if (IC==1){

  
  fraction (f, inside_or_out(y, x, x_bar, x_R, menisc_max, under_x, over_x,
     menisc_x,  under_y,  over_y,  menisc_y,
   under_N,   over_N,  menisc_N));

  refine((inside_or_out_grid_refinement(y, x, x_bar, x_R, menisc_max, under_x, over_x,
    menisc_x,  under_y,  over_y,  menisc_y,
  under_N,   over_N,  menisc_N)) && level < 12);
  }
  free(under_x);
  free(under_y);
  free(over_x);
  free(over_y);
  free(menisc_x);
  free(menisc_y);

  // Strong refinement around the interfacial regions
  //refine (level < maxLevel);
  /*  
  A = sq(x - (filmHeight + 1.0 + 0.5)) + sq(y) < sq(1.0*1.05)  ## Ball at 1.5 above film with radius 1.05
  B = sq(x - (filmHeight + 1.0 + 0.5)) + sq(y) > sq(1.0*0.95)) ## outside Ball at 1.5 above film with radius 0.95
  C  = fabs(x - filmHeight) <= 0.005) ## interface
  && intersect
  */
  // Create active liquid phase as union between drop and film
  if (IC==0){
    refine (((sq(x - (filmHeight - 1.0 + 0.1)) + sq(y) < sq(1.0*1.05) && sq(x - (filmHeight - 1.0 + 0.1)) + sq(y) > sq(1.0*0.95)) || fabs(x - filmHeight) <= 0.005) && level < maxLevel);

    fraction (f, (-sq(1.0) + sq(x - (filmHeight-1.0+0.01)) + sq(y) > 0) && (+ x - filmHeight < 0));
  }
  
  //fraction (f2, (- x + filmHeight));
  
  // Initialise uniform velocity field inside droplet
  printf("after refinement\n");
  foreach()
  {
    
  	u.x[] = 0.0*f[];
        u.y[] = 0.0;
        p[] = 0.0;
	omega[] = 0.0;
  }
}

event adapt (i=5;i++) {

  // Refine only with respect to interfacial shape(s) location and velocity component magnitude
  //printf("late");
  adapt_wavelet ((scalar *){f, u}, (double[]){1e-5, 1e-3, 1e-3}, maxLevel, minLevel);

}
event adapt_early (i=0;i++;i<5) {
  //printf("early");
  // Refine only with respect to interfacial shape(s) location and velocity component magnitude
  adapt_wavelet ((scalar *){f, u}, (double[]){1e-5, 1e-3, 1e-3}, 14, minLevel);

}

/*  event eject_velocity (i++) {
  double jet_vel;
  double temp = max_height;
  double gap = 1e-6;
  
  if (interpolate(f, max_height,0)>0.5){

    for (double j=max_height;j<max_height+0.5; j+= gap){
      if (interpolate(f, j,0)<0.5){
        temp = j-gap;
        break;
      }
    } // int j
  }
  else{
    for (double j=max_height;j>max_height-0.5; j-= gap){
      if (interpolate(f, j,0)>0.5){
        temp = j+gap;
        break;
      }
    } // int j
  }
  
    
  jet_vel = (temp - max_height)/dt;
  max_height = temp;
  
  fprintf(fp_jet_vel,"%g %g %g %g\n", temp,interpolate(u.x,temp,0),jet_vel,t);


}  */

/* event gfsview (t = 0.0; t += 0.1; t <= tEnd) {
    char name_gfs[200];
    sprintf(name_gfs,"Slices/DropImpact-%0.1f.gfs",t);

    FILE* fp_gfs = fopen (name_gfs, "w");
    output_gfs(fp_gfs);
    fclose(fp_gfs);
} */

event saveInterfaces (t = 0.0; t += 0.01; t <= tEnd) {
    char nameInterfaces1[200];

    sprintf(nameInterfaces1,"Interfaces/interfaceDrop-%0.2f.dat",t);

    FILE * fp1 = fopen(nameInterfaces1, "w");
    output_facets (f, fp1);	
    fclose(fp1);

    /* char nameInterfaces2[200];

    sprintf(nameInterfaces2,"Interfaces/interfacePool-%0.1f.dat",t);

    FILE * fp2 = fopen(nameInterfaces2, "w");
    output_facets (f2, fp2);	
    fclose(fp2); */
} 
event loginterface (i += 1) {

  scalar posX[],posY[];
  position (f, posX, {1,0}); // (1,0) indicates the unit vector in the x-direction
  position (f, posY, {0,1}); // (0,1) indicates the unit vector in the y-direction

  fprintf(fp_interface, "%i %g %1.4f %1.4f %1.4f %1.4f %1.4f %1.4f\n", i, t, statsf(f).sum, statsf(posX).min, statsf(posX).max, statsf(posY).min, statsf(posY).max,interpolate(u.x,statsf(posX).max,0));
  fflush(fp_interface);
}

/* event extractPressureData (t = 0.0; t += 0.1) {

  char namePressureData[200];

  sprintf(namePressureData,"Interfaces/customPressureData-%0.1f.dat",t);

  FILE *fpPressureData = fopen(namePressureData, "w");
  
  // Extract pressure data in tailored contact region (for visualisation and analysis inside gas film)
  for(double x = -1.; x < 0.5; x += 0.001){
	for(double y = 0.; y < 1.5; y += 0.001){
	    fprintf (fpPressureData, "%g\n", interpolate (p, x, y));
	}
  }

  fclose(fpPressureData);

} */

// Fluid volume metrics
event droplets (t += 0.01)
{
  scalar m[];
  foreach()
    m[] = f[] > 1e-2;
  int n = tag (m);

  double v[n];
  coord b[n];
  for (int j = 0; j < n; j++)
    v[j] = b[j].x = b[j].y = b[j].z = 0.;
  foreach_leaf()
    if (m[] > 0) {
      int j = m[] - 1;
      v[j] += dv()*f[];
      coord p = {x,y,z};
      foreach_dimension()
	b[j].x += dv()*f[]*p.x;
    }

  #if _MPI
    MPI_Allreduce (MPI_IN_PLACE, v, n, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce (MPI_IN_PLACE, b, 3*n, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  #endif
  for (int j = 0; j < n; j++)
    fprintf (fp_droplets, "%d %g %d %g %g %g\n", i, t,
	     j, v[j], b[j].x/v[j], b[j].y/v[j]);
  fflush (fp_droplets);
}


// Output animations
// event movies (t += 0.01; t <= tEnd){
//
//  char timestring[100];
//  
//  foreach(){
//	omega[] = (u.y[1,0] - u.y[-1,0])/(2.*Delta) - (u.x[0,1] - u.x[0,-1])/(2.*Delta);
//        velnorm[] = sqrt(sq(u.x[]) + sq(u.y[]));
//  	viewingfield[] = 1.0 - f[];// - 0.5*f2[];
//  	mylevel[] = level;
//  }
//
//  view(width=1900, height=1050, fov=7.0, ty = 0.0, quat = { 0, 0, -0.707, 0.707 });
//	
//  clear();
//  draw_vof("f", lw=2);
//  //draw_vof("f2", lw=2);
//  squares("viewingfield", map = cool_warm, min = -0.5, max = 2.5);
//  mirror({0,1}) {
//	draw_vof("f", lw=2);
//	//draw_vof("f2", lw=2);		
//	cells(lw=0.5);
//	squares("mylevel", map = cool_warm, min = minLevel, max = maxLevel);
//  } 
//
//  sprintf(timestring, "t=%2.03f",t);
//  draw_string(timestring, pos=1, lc= { 0, 0, 0 }, lw=2);
//  
//  save ("Animations/ImpactSummary.mp4"); 
//
//  
//  view(width=1900, height=1050, fov=7.0, ty = 0.0, quat = { 0, 0, -0.707, 0.707 });;
//  clear();
//  
//  draw_vof("f", lw=2);
//  draw_vof("f2", lw=2);
//  squares("u.x", map = cool_warm, min = -1., max = 0.5,);
//  mirror({0,1}) {
//	draw_vof("f", lw=2);
//	draw_vof("f2", lw=2);		
//	squares("u.y", map = cool_warm, min = -0.5, max = 2.);
//  } 
//
//  sprintf(timestring, "t=%2.03f",t);
//  draw_string(timestring, pos=1, lc= { 0, 0, 0 }, lw=2);
//  
//  save ("Animations/ImpactVelocities.mp4"); 
//
//  view(width=1900, height=1050, fov=7.0, ty = 0.0, quat = { 0, 0, -0.707, 0.707 });
//  clear();
//  
//  draw_vof("f", lw=2);
//  draw_vof("f2", lw=2);
//  squares("omega", map = cool_warm, min = -3., max = 3.);
//  mirror({0,1}) {
//	draw_vof("f", lw=2);
//	draw_vof("f2", lw=2);	
//	squares("p", map = cool_warm, min = -0.25, max = 4.0);
//  } 
//
//  sprintf(timestring, "t=%2.03f",t);
//  draw_string(timestring, pos=1, lc= { 0, 0, 0 }, lw=2);
//  
//  save ("Animations/ImpactPVort.mp4"); 
//  
//} 

event logstats (t += 0.01; t <= tEnd) {

    timing s = timer_timing (perf.gt, i, perf.tnc, NULL);
 
    // Output i, timestep, no of cells, real time elapsed, cpu time
    fprintf(fp_stats, "i: %i t: %g dt: %g #Cells: %ld Wall clock time (s): %g CPU time (s): %g \n", i, t, dt, grid->n, perf.t, s.cpu);
    fflush(fp_stats);
}

/*###########################################################################*/
/*###########################################################################*/
/*###########################################################################*/
/*###########################################################################*/
/*###########################################################################*/
/*###########################################################################*/
/*###########################################################################*/

/******************************************************************************/

void spline_linear_val ( int ndata, double tdata[], double ydata[], 
  double tval, double *yval, double *ypval )

/******************************************************************************/
/*
  Purpose:

    SPLINE_LINEAR_VAL evaluates a piecewise linear spline at a point.

  Discussion:

    Because of the extremely simple form of the linear spline,
    the raw data points ( TDATA(I), YDATA(I)) can be used directly to
    evaluate the spline at any point.  No processing of the data
    is required.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    24 February 2004

  Author:

    John Burkardt

  Parameters:

    Input, int NDATA, the number of data points defining the spline.

    Input, double TDATA[NDATA], YDATA[NDATA], the values of the independent
    and dependent variables at the data points.  The values of TDATA should
    be distinct and increasing.

    Input, double TVAL, the point at which the spline is to be evaluated.

    Output, double *YVAL, *YPVAL, the value of the spline and its first
    derivative dYdT at TVAL.  YPVAL is not reliable if TVAL is exactly
    equal to TDATA(I) for some I.
*/
{
  int left;
  int right;
/*
  Find the interval [ TDATA(LEFT), TDATA(RIGHT) ] that contains, or is
  nearest to, TVAL.
*/
  r8vec_bracket ( ndata, tdata, tval, &left, &right );
/*
  Now evaluate the piecewise linear function.
*/
  *ypval = ( ydata[right-1] - ydata[left-1] ) 
         / ( tdata[right-1] - tdata[left-1] );

  *yval = ydata[left-1] +  ( tval - tdata[left-1] ) * (*ypval);

  return;
}

/******************************************************************************/

void r8vec_bracket ( int n, double x[], double xval, int *left,
  int *right )

/******************************************************************************/
/*
  Purpose:

    R8VEC_BRACKET searches a sorted array for successive brackets of a value.

  Discussion:

    An R8VEC is a vector of R8's.

    If the values in the vector are thought of as defining intervals
    on the real line, then this routine searches for the interval
    nearest to or containing the given value.

    It is always true that RIGHT = LEFT+1.

    If XVAL < X[0], then LEFT = 1, RIGHT = 2, and
      XVAL   < X[0] < X[1];
    If X(1) <= XVAL < X[N-1], then
      X[LEFT-1] <= XVAL < X[RIGHT-1];
    If X[N-1] <= XVAL, then LEFT = N-1, RIGHT = N, and
      X[LEFT-1] <= X[RIGHT-1] <= XVAL.

    For consistency, this routine computes indices RIGHT and LEFT
    that are 1-based, although it would be more natural in C and
    C++ to use 0-based values.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    31 May 2009

  Author:

    John Burkardt

  Parameters:

    Input, int N, length of input array.

    Input, double X[N], an array that has been sorted into ascending order.

    Input, double XVAL, a value to be bracketed.

    Output, int *LEFT, *RIGHT, the results of the search.
*/
{
  int i;

  for ( i = 2; i <= n - 1; i++ )
  {
    if ( xval < x[i-1] )
    {
      *left = i - 1;
      *right = i;
      return;
    }

   }

  *left = n - 1;
  *right = n;

  return;
}