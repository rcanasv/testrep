#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define NAME_LENGTH 256
#define SKIP    fread (&dummyi,  sizeof(int), 1, f);  
#define PSKIP   printf ("%d\n", dummyi);

struct pdata
{
  double   Pos[3];
  double   Vel[3];
  double   Mass;
  int      Id;
  int      Type;
  double   Age;
  double   Metal;
  double   Chem;
  double   Radius;
};

struct hmProps
{
  int      nparts;
  int   *  partIDs;
  int      id;
  int      timestep;
  int      level;
  int      parentid;
  int      idfirstchild;
  int      nsubs;
  int      idnextgal;
  float    mass;
  float    pos[3];
  float    vel[3];
  float    ang[3];
  float    rsize;
  float    eigval[3];
  float    Ekin;
  float    Epot;
  float    Etot;
  float    lambda;
  float    vdisp;
  float    bvdisp;
  float    bmass;
  float    rvir;
  float    mvir;
  float    tvir;
  float    cs;
  float    nfw_rho;
  float    nfw_rad;
  int      nbins;
  float *  rbin;
  float *  sbin;
  float    mhalf;
  float    rhalf;
};

struct hmOutput
{
  char              tbricks     [NAME_LENGTH];
  char              galpath     [NAME_LENGTH];
  char              snappath    [NAME_LENGTH];
  char              snapprefix  [NAME_LENGTH];
  int               nparts;
  float             massres;
  float             aexp;
  float             OmegaM;
  float             AgeUniv;
  int               nsubs;
  int               nstructs;
  int               ntotal;
  struct hmProps  * strctProps;
  struct pdata   ** strctParts;
};


void read_halomaker_treebricks (struct hmOutput * hmo, int iReadParticles, int iGetProfile);
void read_gal_file (char * path, struct hmProps * strct, struct pdata ** p);

int rad_compare (const void * a, const void * b);

double friedman(double Omega0, double OmegaL, double OmegaK, double alpha, double axp_min, double ** axp_out, double ** hexp_out, double ** tau_out, double ** t_out, int ntable);
double dadtau  (double axp_tau, double Omega0, double OmegaL, double OmegaK);
double dadt    (double axp_t, double Omega0, double OmegaL, double OmegaK);



int main (int argc, char ** argv)
{
  int i, j, k;
  FILE * f;
  char fname [NAME_LENGTH];
  
  struct hmOutput hm;
  
  strcpy(hm.tbricks,     argv[1]);
  strcpy(hm.galpath,     argv[2]);
  strcpy(hm.snappath,    argv[3]);
  strcpy(hm.snapprefix,  argv[4]);
  
  read_halomaker_treebricks (&hm, 0, 0);
  
  
  
  hm.strctParts = malloc (hm.ntotal*sizeof(struct pdata *));
  for (i = 1; i <= hm.ntotal; i++)
  {    
    read_gal_file (hm.galpath, &hm.strctProps[i], &hm.strctParts[i]);
    for (j = 0; j < hm.strctProps[i].nparts; j++)
    {
      hm.strctParts[i][j].Pos[0] *= 1000;
      hm.strctParts[i][j].Pos[1] *= 1000;
      hm.strctParts[i][j].Pos[2] *= 1000;
      hm.strctParts[i][j].Mass   *= 1e11;
    }
  }
  
  
//   //
//   //   move_gal_to_cm ()
//   //   
//   
//   double cm[3];
//   double mtot;
//   double rhalf;
//   double mhalf;
//   
//   for (i = 1; i <= hm.ntotal; i++)
//   {
//     // First calculate center of mass
//     cm[0] = 0.0;
//     cm[1] = 0.0;
//     cm[2] = 0.0;
//     mtot  = 0.0;
//     
//     for (j = 0; j < hm.strctProps[i].nparts; j++)
//     {
//       cm[0] += hm.strctParts[i][j].Pos[0] * hm.strctParts[i][j].Mass;
//       cm[1] += hm.strctParts[i][j].Pos[1] * hm.strctParts[i][j].Mass;
//       cm[2] += hm.strctParts[i][j].Pos[2] * hm.strctParts[i][j].Mass;
//       mtot  += hm.strctParts[i][j].Mass;
//     }
//     cm[0] /= mtot;
//     cm[1] /= mtot;
//     cm[2] /= mtot;
//     
//     // Shift particles to center of Mass
//     for (j = 0; j < hm.strctProps[i].nparts; j++)
//     {
//       hm.strctParts[i][j].Pos[0] -= cm[0];
//       hm.strctParts[i][j].Pos[1] -= cm[1];
//       hm.strctParts[i][j].Pos[2] -= cm[2];
//       hm.strctParts[i][j].Radius = hm.strctParts[i][j].Pos[0]*hm.strctParts[i][j].Pos[0] + \
//                                    hm.strctParts[i][j].Pos[1]*hm.strctParts[i][j].Pos[1] + \
//                                    hm.strctParts[i][j].Pos[2]*hm.strctParts[i][j].Pos[2];
//     }
//     
//     // Sort by radius
//     qsort (hm.strctParts[i], hm.strctProps[i].nparts, sizeof(struct pdata), rad_compare);
//     
//     // Calculate half mass radius
//     hm.strctProps[i].rhalf = 0;
//     hm.strctProps[i].mhalf = 0;
//     
//     for (j = 0; j < hm.strctProps[i].nparts; j++)
//     {
//       hm.strctProps[i].mhalf += hm.strctParts[i][j].Mass;
//       if (hm.strctProps[i].mhalf > (0.5*mtot))
//       {
//         hm.strctProps[i].rhalf = sqrt(hm.strctParts[i][j].Radius);
//         hm.strctProps[i].mass  = mtot;
//         break;
//       }
//     }
//   }
// 
//  
//   //
//   // Get Mass Size - Relation
//   //
//   
//   // Write mass size relation file
//   sprintf (fname, "MassSize_hmkr_%03d.dat", hm.strctProps[0].timestep);
//   f = fopen (fname, "w");
//   for (i = 1; i <= hm.ntotal; i++)
//     fprintf (f, "%e   %e   %e\n",hm.strctProps[i].rhalf, hm.strctProps[i].mhalf, hm.strctProps[i].mass);
//   fclose (f);
//   
// 

  
  
//   //
//   //  Compute SFR
//   //
//   double boxlen;
//   double time;
//   double aexp;
//   double H0;
//   double Omega_m;
//   double Omega_l;
//   double Omega_b;
//   char name [NAME_LENGTH];
//   char dummys [NAME_LENGTH];
//   double t;
//   char   buffer       [NAME_LENGTH];
// 
//   sprintf (name, "%s/info_%05d.txt", hm.snappath, hm.strctProps[1].timestep);  
//   FILE * ff = fopen(name,"r");
//   for (i = 0; i < 7; i++)
//     fgets(buffer, 100, ff);
//   fgets(buffer, 100, ff);   sscanf (buffer, "%s %s %lf", dummys, dummys, &boxlen);
//   fgets(buffer, 100, ff);   sscanf (buffer, "%s %s %lf", dummys, dummys, &time);
//   fgets(buffer, 100, ff);   sscanf (buffer, "%s %s %lf", dummys, dummys, &aexp);
//   fgets(buffer, 100, ff);   sscanf (buffer, "%s %s %lf", dummys, dummys, &H0);
//   fgets(buffer, 100, ff);   sscanf (buffer, "%s %s %lf", dummys, dummys, &Omega_m);
//   fgets(buffer, 100, ff);   sscanf (buffer, "%s %s %lf", dummys, dummys, &Omega_l);
//   fclose (ff);
// 
//   double * axp_frw;
//   double * hexp_frw;
//   double * tau_frw;
//   double * t_frw;
//   int      n_frw = 1000;
//   double   time_tot = friedman(Omega_m, Omega_l, 0.0, 1e-6, 1e-3, &axp_frw, &hexp_frw, &tau_frw, &t_frw, n_frw);
// 
//   // Find neighbouring conformal time
//   i = 1;
//   while (tau_frw[i] > time  && i < n_frw)
//     i = i+1;
//   
//   // Interpolate time
//   double time_simu = t_frw[i]   * (time - tau_frw[i-1]) / (tau_frw[i]   - tau_frw[i-1]) + \
//                     t_frw[i-1] * (time - tau_frw[i])   / (tau_frw[i-1] - tau_frw[i]);
//                   
//   printf ("Time simu    %lf\n", (time_tot + time_simu) / (H0*1e5/3.08e24) / (365*24*3600*1e9));
//   printf ("Hubble time  %lf\n", time_tot / (H0*1e5/3.08e24) / (365*24*3600*1e9));
//   
//   printf ("i               %d\n", i);
//   printf ("time            %e\n", time);
//   printf ("time_tot        %e\n", time_tot);
//   printf ("time_simu       %e\n", time_simu);
//   printf ("t_tot + t_simu  %e\n", time_tot + time_simu);
// 
//   for (i = 1; i <= hm.ntotal; i++)
//   {
//     for (j = 0; j < hm.strctProps[i].nparts; j++)
//     {
//       k = 1;
//       while (tau_frw[k] > hm.strctParts[i][j].Age  && k < n_frw)
//         k++;
//       
//       t = t_frw[k]   * (hm.strctParts[i][j].Age - tau_frw[k-1]) / (tau_frw[k]   - tau_frw[k-1]) + \
//           t_frw[k-1] * (hm.strctParts[i][j].Age - tau_frw[k])   / (tau_frw[k-1] - tau_frw[k]);
//       hm.strctParts[i][j].Age = (time_simu - t) / (H0*1e5/3.08e24) / (365*24*3600.0);
//     }
//   }
//   free (axp_frw);
//   free (hexp_frw);
//   free (tau_frw);
//   free (t_frw);
  
//   char bfr [NAME_LENGTH];
//   sprintf (bfr, "sSFR_hmkr_%d.dat", hm.strctProps[1].timestep);
//   FILE * fgal = fopen (bfr, "w");    
//   double Mnew;
//   int nage;
//   double Agesum;
//   double mtot;
//   double Dt = 100e6;
//   for (i = 1; i <= hm.ntotal; i++)
//   {
//     Mnew = 0;
//     mtot = 0;
//     nage = 0;
//     Agesum = 0;
//     for (j = 0; j < hm.strctProps[i].nparts; j++)
//     {
//       if (hm.strctParts[i][j].Age < Dt && hm.strctParts[i][j].Age > 0.0)
// //           if (P[i][j].Age < Dt)
//       {
//         Mnew += hm.strctParts[i][j].Mass;
//         nage++;
//         Agesum += hm.strctParts[i][j].Age;
//       }
//       mtot += hm.strctParts[i][j].Mass;
//     }
//     if (nage)
//       fprintf (fgal, "%d   %e   %e   %e   %e\n", i, hm.strctProps[i].mass, mtot, Mnew/Dt, Agesum/(double)nage);
//   }
//   fclose (fgal);
  
  
  
  
//   for (i = 1; i <= hm.ntotal; i++)
//     free (hm.strctParts[i]);
//   free (hm.strctProps);
  
   
//   free_hmOutput (&hm);
    
  return 0;
}


int rad_compare (const void * a, const void * b)
{
  struct pdata * Part1 = (struct pdata *)a;
  struct pdata * Part2 = (struct pdata *)b;
  
  if (Part1->Radius > Part2->Radius)
    return 1;
  if (Part1->Radius == Part2->Radius)
    return 0;
  if (Part1->Radius < Part2->Radius)
    return -1;  
}


//
//  Read HaloMaker TreeBricks
//
void read_halomaker_treebricks (struct hmOutput * hmo, int iReadParticles, int iGetProfile)
{
  int i, j, k;
  
  int    dummyi;
  float  dummyf;
  
  FILE * f;
  
  f = fopen (hmo->tbricks, "r");
  
  //
  // Read HaloMaker tree_bricks header
  //
  SKIP  fread (&hmo->nparts,   sizeof(int),   1, f);   SKIP
  SKIP  fread (&hmo->massres,  sizeof(float), 1, f);   SKIP
  SKIP  fread (&hmo->aexp,     sizeof(float), 1, f);   SKIP
  SKIP  fread (&hmo->OmegaM,   sizeof(float), 1, f);   SKIP
  SKIP  fread (&hmo->AgeUniv,  sizeof(float), 1, f);   SKIP
  SKIP  fread (&hmo->nstructs, sizeof(int),   1, f);  
        fread (&hmo->nsubs,    sizeof(int),   1, f);   SKIP

  hmo->ntotal = hmo->nstructs + hmo->nsubs;
  
  hmo->strctProps = (struct hmProps *) malloc ((hmo->ntotal+1) * sizeof(struct hmProps));
  
  //
  // Loop over structures
  //
  for (i = 1; i <= hmo->ntotal; i++)
  {
    // Number of particles
    SKIP  
    fread (&hmo->strctProps[i].nparts, sizeof(int), 1, f);
    SKIP
    
    // List of particles in galaxy i
    SKIP
    if (iReadParticles)
    {
      hmo->strctProps[i].partIDs = (int *) malloc (hmo->strctProps[i].nparts * sizeof(int));
      fread (hmo->strctProps[i].partIDs,  sizeof(int), hmo->strctProps[i].nparts, f);
    }
    else
      fseek (f, dummyi, SEEK_CUR);
    SKIP 
    
    // Galaxy Identity
    SKIP
    fread (&hmo->strctProps[i].id, sizeof(int), 1, f);
    SKIP
    
    // TimeStep Number
    SKIP
    fread (&hmo->strctProps[i].timestep, sizeof(int), 1, f);
    SKIP

    // Level, ParentID, ID1stchild, nsubs, idnextgal
    SKIP
    fread (&hmo->strctProps[i].level, sizeof(int), 5, f);
    SKIP
    
    // Galaxy Mass
    SKIP
    fread (&hmo->strctProps[i].mass, sizeof(float), 1, f);
    SKIP
    
    // Galaxy Position
    SKIP
    fread (&hmo->strctProps[i].pos[0], sizeof(float), 3, f);
    SKIP
    
    // Galaxy Velocity
    SKIP
    fread (&hmo->strctProps[i].vel[0], sizeof(float), 3, f);
    SKIP
    
    // Galaxy Angular Momentum
    SKIP
    fread (&hmo->strctProps[i].ang[0], sizeof(float), 3, f);
    SKIP
    
    // Distance to most distant particle, inertia tensor eigvals
    SKIP
    fread (&hmo->strctProps[i].rsize, sizeof(float), 4, f);
    SKIP
    
    // Kinetic Potential and Total Energy
    SKIP
    fread (&hmo->strctProps[i].Ekin, sizeof(float), 3, f);
    SKIP
    
    // Spin Parameter
    SKIP
    fread (&hmo->strctProps[i].lambda, sizeof(float), 1, f);
    SKIP
    
    // Velocity dispersion, bulge velocity disp, bulge mass
    SKIP
    fread (&hmo->strctProps[i].vdisp, sizeof(float), 3, f);
    SKIP
    
    // Virial Radius, Virial Mass, Virial Temp, Sound Speed
    SKIP
    fread (&hmo->strctProps[i].rvir, sizeof(float), 4, f);
    SKIP
    
    // Central Density (NFW), Characteristic Radius (NFW)
    SKIP
    fread (&hmo->strctProps[i].nfw_rho, sizeof(int), 2, f);
    SKIP
    
    // Stellar surface density profiles
    SKIP
    fread (&hmo->strctProps[i].nbins, sizeof(int), 1, f);
    SKIP
    
    SKIP
    if (iGetProfile)
    {
      hmo->strctProps[i].rbin = (float *) malloc (hmo->strctProps[i].nbins * sizeof(float));
      fread (hmo->strctProps[i].rbin, sizeof(float), hmo->strctProps[i].nbins, f);
    }
    else
      fseek (f, dummyi, SEEK_CUR);
    SKIP
      
    SKIP
    if (iGetProfile)
    {
      hmo->strctProps[i].sbin = (float *) malloc (hmo->strctProps[i].nbins * sizeof(float));
      fread (hmo->strctProps[i].sbin, sizeof(float), hmo->strctProps[i].nbins, f);
    }
    else
      fseek (f, dummyi, SEEK_CUR);
    SKIP
  }
  fclose (f);
}

//
// Read Galfile
//
void read_gal_file (char * path, struct hmProps * strct, struct pdata ** p)
{  
  int     dummyi;
  int     i, j;
  char    filename [NAME_LENGTH];
  FILE *  f;
  
  
  int     gal_number;
  int     gal_level;
  double  gal_mass;
  double  gal_pos[3];
  double  gal_vel[3];
  double  gal_ang[3];
  int     nlist;

  sprintf (filename, "%s/gal_stars_%07d", path, strct->id);
  if ((f = fopen (filename, "r")) == NULL)
  {
    printf ("Can't open file named   %s \n", filename);
    exit(0);
  }

  SKIP  fread (&gal_number, sizeof(int),    1, f);  SKIP
  SKIP  fread (&gal_level,  sizeof(int),    1, f);  SKIP  
  SKIP  fread (&gal_mass,   sizeof(double), 1, f);  SKIP
  SKIP  fread (&gal_pos,    sizeof(double), 3, f);  SKIP
  SKIP  fread (&gal_vel,    sizeof(double), 3, f);  SKIP
  SKIP  fread (&gal_ang,    sizeof(double), 3, f);  SKIP
  SKIP  fread (&nlist,      sizeof(int),    1, f);  SKIP

  if ((strct->id     != gal_number) ||  \
      (strct->nparts != nlist)          \
     )
  {
     printf ("Treebricks and GAL_ Properties don't match\n");
     printf ("Exiting...\n");
     exit(0);
  }
  
  if (!(p[0] = malloc (nlist * sizeof(struct pdata))))
  {
    printf ("Cannot allocate memory for particle information\n");
    printf ("Exiting\n");
    exit (0);
  }

  for (j = 0; j < 3; j++)
  {
    SKIP
    if (dummyi == 8 * nlist)
      for (i = 0; i < nlist; i++)
        fread (&p[0][i].Pos[j], sizeof(double), 1, f);
    SKIP
  }
  
  for (j = 0; j < 3; j++)
  {
    SKIP
    if (dummyi == 8 * nlist)
      for (i = 0; i < nlist; i++)
        fread (&p[0][i].Vel[j], sizeof(double), 1, f);
    SKIP
  }        

  SKIP
  if (dummyi == 8 * nlist)
    for (i = 0; i < nlist; i++)
      fread (&p[0][i].Mass, sizeof(double), 1, f);
  SKIP
  
  SKIP
  if (dummyi == 4 * nlist)
    for (i = 0; i < nlist; i++)
      fread (&p[0][i].Id, sizeof(int), 1, f);
  SKIP
  
  SKIP
  if (dummyi == 8 * nlist)
    for (i = 0; i < nlist; i++)
      fread (&p[0][i].Age, sizeof(double), 1, f);
  SKIP

  SKIP
  fseek (f, dummyi, SEEK_CUR);
//   if (dummy == 8 * nlist)
//     for (i = 0; i < nlist; i++)
//       fread (&p[0][i].Metal, sizeof(double), 1, f);
  SKIP
  
  SKIP
  fseek (f, dummyi, SEEK_CUR);
//   if (dummy == 8 * nlist)
//     for (i = 0; i < nlist; i++)
//       fread (&p[0][i].Chem, sizeof(double), 1, f);
  SKIP
  
  fclose(f); 
}


//
//  Friedman Lookup tables
//
double friedman(double Omega0, double OmegaL, double OmegaK, double alpha, double axp_min, double ** axp_out, double ** hexp_out, double ** tau_out, double ** t_out, int ntable)
{
  /*
  ! ######################################################!
  ! This subroutine assumes that axp = 1 at z = 0 (today) !
  ! and that t and tau = 0 at z = 0 (today).              !
  ! axp is the expansion factor, hexp the Hubble constant !
  ! defined as hexp=1/axp*daxp/dtau, tau the conformal    !
  ! time, and t the look-back time, both in unit of 1/H0. !
  ! alpha is the required accuracy and axp_min is the     !
  ! starting expansion factor of the look-up table.       !
  ! ntable is the required size of the look-up table.     !
  ! ######################################################!
  */
  double axp_tau = 1.0;
  double axp_t   = 1.0;
  double tau     = 0.0;
  double t       = 0.0;
  double age_tot;
  
  double dtau;
  double dt;
  double axp_tau_pre;
  double axp_t_pre;
  
  int    nstep   = 0;
  int nskip;
  int nout;
  
  *axp_out  = (double *) malloc (ntable * sizeof(double));
  *hexp_out = (double *) malloc (ntable * sizeof(double));
  *tau_out  = (double *) malloc (ntable * sizeof(double));
  *t_out    = (double *) malloc (ntable * sizeof(double));
  
  while ((axp_tau >= axp_min) || (axp_t >= axp_min))
  {
    nstep++;
    dtau = alpha * axp_tau / dadtau(axp_tau, Omega0, OmegaL, OmegaK);
    axp_tau_pre = axp_tau - dadtau(axp_tau, Omega0, OmegaL, OmegaK) * dtau / 2.0;
    axp_tau = axp_tau - dadtau(axp_tau_pre, Omega0, OmegaL, OmegaK) * dtau;
    tau = tau - dtau;
    
    dt = alpha * axp_t / dadt(axp_t, Omega0, OmegaL, OmegaK);
    axp_t_pre = axp_t - dadt(axp_t, Omega0, OmegaL, OmegaK) * dt / 2.0;
    axp_t = axp_t - dadt(axp_t_pre, Omega0, OmegaL, OmegaK) * dt;
    t = t - dt;
  }
  
  age_tot =-t;
//   printf ("Age of the universe (in unit of 1/H0)=%e \n", -t);
  
  nskip = nstep / ntable;
  
  axp_t   = 1.0;
  axp_tau = 1.0;
  tau     = 0.0;
  t       = 0.0;
  
  nstep = 0;
  nout  = 0;
  
  t_out   [0][nout] = t;
  tau_out [0][nout] = tau;
  axp_out [0][nout] = axp_tau;
  hexp_out[0][nout] = dadtau (axp_tau, Omega0, OmegaL, OmegaK) / axp_tau;

  
  while ((axp_tau >= axp_min) || (axp_t >= axp_min))
  {
    nstep++;
    dtau = alpha * axp_tau / dadtau (axp_tau, Omega0, OmegaL, OmegaK);
    axp_tau_pre = axp_tau - dadtau(axp_tau, Omega0, OmegaL, OmegaK) * dtau/2.0;
    axp_tau = axp_tau - dadtau(axp_tau_pre, Omega0, OmegaL, OmegaK) * dtau;
    tau = tau - dtau;
    
    dt = alpha * axp_t / dadt(axp_t, Omega0, OmegaL, OmegaK);
    axp_t_pre = axp_t - dadt(axp_t, Omega0, OmegaL, OmegaK) * dt / 2.0;
    axp_t = axp_t - dadt(axp_t_pre, Omega0, OmegaL, OmegaK) * dt;
    t = t -dt;
    
    if ((nstep%nskip) == 0)
    {
      nout = nout + 1;
      t_out   [0][nout] = t;
      tau_out [0][nout] = tau;
      axp_out [0][nout] = axp_tau;
      hexp_out[0][nout] = dadtau(axp_tau, Omega0, OmegaL, OmegaK) / axp_tau;
    }
  }
  
  t_out   [0][ntable-1] = t;
  tau_out [0][ntable-1] = tau;
  axp_out [0][ntable-1] = axp_tau;
  hexp_out[0][ntable-1] = dadtau(axp_tau, Omega0, OmegaL, OmegaK) / axp_tau;
  
  return age_tot;
}

//
//  DADTAU
//
double dadtau (double axp_tau, double Omega0, double OmegaL, double OmegaK)
{
  return sqrt(axp_tau*axp_tau*axp_tau * (Omega0 + OmegaL*axp_tau*axp_tau*axp_tau + OmegaK*axp_tau));
}

//
//  DADT
//
double dadt (double axp_t, double Omega0, double OmegaL, double OmegaK)
{
  return sqrt((1.0/axp_t) * (Omega0 + OmegaL*axp_t*axp_t*axp_t + OmegaK*axp_t));
}


