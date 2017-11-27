#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
// #include <hdf5.h>


#define NAME_LENGTH 256
#define LONG_LENGTH 3000
//#define NUMFILES 41

int     ndim;
int     gal_number;
int     gal_level;
double  gal_mass;
double  gal_pos[3];
double  gal_vel[3];
double  gal_ang[3];
int     nlist;

int * extended_oIndex;
int * extended_IdStruct;
int * extended_IdHost;
int * extended_IdIGM;

double ** ramses_pos;
double ** ramses_vel;
double ** ramses_met;
double  * ramses_age;
double  * ramses_mass;
int     * ramses_id;
int     * ramses_lvl;


struct pdata
{
  double   Pos[3];
  double   Vel[3];
  double   Mass;
  int      Id;
  double   Age;
  double   Metal;
  double   Chem;
  int      Type;
};

struct io_header_1
{
  int      npart[6];
  double   mass[6];
  double   time;
  double   redshift;
  int      flag_sfr;
  int      flag_feedback;
  int      npartTotal[6];
  int      flag_cooling;
  int      num_files;
  double   BoxSize;
  double   Omega0;
  double   OmegaLambda;
  double   HubbleParam;
  char     fill[256 - 6*4 - 6*8 - 2*8 - 2*4 - 6*4 - 2*4 - 4*8];      /* fills to 256 Bytes */
} header1;

struct pgadget
{
  float    Pos[3];
  float    Vel[3];
  float    Mass;
  int      Id;
  int      Type;
  float    Age;
  float    Metal;
  float    Chem;
} ** P, * p;

int * ID;
int NumPart;
int outType;

char   buffer       [NAME_LENGTH];
char   longbuffer   [LONG_LENGTH];
char   stfprefix    [NAME_LENGTH];
char   galprefix    [NAME_LENGTH];
char   outprefix    [NAME_LENGTH];
char   directory    [NAME_LENGTH];
char   format       [NAME_LENGTH];
char   output_fname [NAME_LENGTH];


struct objProps
{
  int      ID;
  int      DirectHostID;
  int      HostID;
  int      NumSubs;
  int      Type;
  int      NumPart;
  double   TotMass;
  double   Pos[3];
  double   Vel[3];
  double   Efrac;
  double   Rsize;
  double   RHalfMass;
  double   Vmax;
  double   Rvmax;
  double   Vdisp;
  double   Lambda;
  double   L[3];
  int    * SubIDs;
  int      NumFiles;
  int    * FilesOfGroup;
  int      NumProg;
  int    * ProgIDs;
  double * ProgMrrts;
  int      dummy;
};


struct stfOutput
{
  char              prefix[NAME_LENGTH];
  int               nstruct;
  int               nprocs;
  int               iprops;
  int               iparts;
  struct objProps * strctProps;
  struct pgadget ** strctParts;
};


int NUMFILES;

int load_stf_extended_output (char * prefix, int filenum);

void read_properties_file (struct stfOutput * tmpstf);

void read_gal_file (char * filename);
void read_gadget_snapshot(char * snapshot);
void read_ramses_snapshot(char * dir, char * snapshot, int numfile);

void write_snapshot(struct pgadget * pg, int particles, struct io_header_1 header, char * output_name);

void free_extended_arrays (void);
void free_ramses_arrays   (void);

int read_stf_filesofgroup (char * prefix, int strct_id, int ** files_of_strct);
int get_n_num_from_string (char * strng, int n_num, int ** nums);

int load_treefrog (char * tffile, int strct_id, int ** prog_ids, float ** prog_mrrts);

void init_stfOutput (struct stfOutput * tmpstf);
void free_stfOutput (struct stfOutput * tmpstf);
void fill_SubIDS    (struct stfOutput * tmpstf);
void fill_ProgIDs (struct stfOutput * tmpstf, char * tffile);

double friedman(double Omega0, double OmegaL, double OmegaK, double alpha, double axp_min, double ** axp_out, double ** hexp_out, double ** tau_out, double ** t_out, int ntable);
double dadtau (double axp_tau, double Omega0, double OmegaL, double OmegaK);
double dadt (double axp_t, double Omega0, double OmegaL, double OmegaK);

int main (int argc, char ** argv)
{
    if (argc < 5)
    {
      printf ("ERROR: Input Parameters missing \n");
      printf ("Usage:    %s   stf_prefix  directory   gal_prefix  output_prefix  format\n", argv[0]);
      exit (0);
    }
    
    FILE * f;
    int    i, j, k;
    
    strcpy (stfprefix, argv[1]);
    strcpy (directory, argv[2]);
    strcpy (galprefix, argv[3]);
    strcpy (outprefix, argv[4]);
    strcpy (format,    argv[5]);
    NUMFILES = atoi (argv[6]);
    int ID = atoi (argv[7]);
    

    if ( (strcmp(format, "gadget") != 0) && (strcmp(format,"ramses") != 0) )
    {
      printf ("Incorrect specified format  %s, exiting...\n", format);
      exit (0);
    }

    struct stfOutput snap1;
    
    init_stfOutput (&snap1);
    
    sprintf (snap1.prefix, "%s", stfprefix);
    
    read_properties_file (&snap1);
    
    
    int * files;
    int nfiles;
    nfiles = read_stf_filesofgroup (stfprefix, ID, &files);
    
    printf ("%d\n", snap1.strctProps[ID].NumPart);
    
    if (!(p = malloc (snap1.strctProps[ID].NumPart * sizeof(struct pgadget))))
    {
      printf ("Cannot allocate memory for particle information for structure %d\n", i);
      printf ("Exiting\n");
      exit (0);
    }
    
    int ninextended = 0;
    int dummystruct = 0;
    int dummyindex = 0;
    int acum = 0;
    
    for (i = 0; i < nfiles; i++)
    {
      printf ("file %d/%d  %d\n", i+1, nfiles, files[i]);
      ninextended = load_stf_extended_output (stfprefix, files[i]);
      printf ("%d\n", ninextended);
      if (ninextended)
      {
        read_ramses_snapshot (directory, galprefix, files[i]);
        for (j = 0; j < ninextended; j++)
        {
          dummystruct = extended_IdStruct[j];
          dummyindex  = extended_oIndex[j];
          if (dummystruct == ID)
          {
            p[acum].Pos[0] = ramses_pos[0][dummyindex];
            p[acum].Pos[1] = ramses_pos[1][dummyindex];
            p[acum].Pos[2] = ramses_pos[2][dummyindex];
            p[acum].Vel[0] = ramses_vel[0][dummyindex];
            p[acum].Vel[1] = ramses_vel[1][dummyindex];
            p[acum].Vel[2] = ramses_vel[2][dummyindex];
            p[acum].Mass   = ramses_mass  [dummyindex];
            p[acum].Id     = ramses_id    [dummyindex];
            p[acum].Age    = ramses_age   [dummyindex];
            p[acum].Type = 4;
            acum++;
          }
        }
        printf ("%d\n", acum);
        free_extended_arrays ();
        free_ramses_arrays ();
      }
    }
    
    write_snapshot(p, acum, header1, outprefix);

    
//     exit(0);


//     int    * prog_ids;
//     float  * prog_mrrts;
//     int      nprogs;
//     char tfname [NAME_LENGTH];
//     sprintf (tfname, "tftest");
//     nprogs = load_treefrog (tfname, 1000001, &prog_ids, &prog_mrrts);
//     exit(0);


    

//     struct stfOutput snap2;
    
//     init_stfOutput (&snap2);
    
//     sprintf (snap1.prefix, "782_gals");
//     sprintf (snap2.prefix, "772_gals");
    
    
    
//     char fname [NAME_LENGTH];
//     sprintf (fname, "MassSize_veloci_%03d.dat", atoi(stfprefix));
//     f = fopen(fname, "w");
//     for (i = 1; i <= snap1.nstruct; i++)
//       if (snap1.strctProps[i].Type > 7)
//       {
//         fprintf (f, "%e  ", snap1.strctProps[i].RHalfMass);
//         fprintf (f, "%e  ", snap1.strctProps[i].TotMass);
//         fprintf (f, "\n");
//       }
//     fclose (f);
    
//     printf ("Properties file read\n");
    

//     read_properties_file (&snap2);
    
//     fill_SubIDS (&snap1);
    
//     fill_SubIDS (&snap2);
    
//     fill_ProgIDs (&snap1, "tftest");
//     
//     int    id1, id2, tmpid;
//     double m1, m2;
//     double ratio;
//     int    nmajor, nminor, nsmooth;
//     double maxratio;
//     double tmmajor, tmminor, tmsmooth;
//     
//     
//     FILE * fmergers = fopen ("mergers.dat", "w");
//     
//     for (i = 1; i <= snap1.nstruct; i++)
//     {
//       if (snap1.strctProps[i].Type != 7)
//       {
//         nmajor  = 0;
//         nminor  = 0;
//         nsmooth = 0;
//         
//         m1       = 0.0;
//         ratio    = 0.0;
//         maxratio = 0.0;
//         
//         tmmajor  = 0.0;
//         tmminor  = 0.0;
//         tmsmooth = 0.0;
//         
//         if (snap1.strctProps[i].NumProg > 1)
//         {
//           for (j = 0; j < snap1.strctProps[i].NumProg; j++)
//           {
//             tmpid = snap1.strctProps[i].ProgIDs[j];
//             if (snap2.strctProps[tmpid].TotMass > m1)
//             {
//               m1  = snap2.strctProps[tmpid].TotMass;
//               id1 = tmpid;
//             }
//           }
//           
//           for (j = 0; j < snap1.strctProps[i].NumProg; j++)
//           {
//             id2 = snap1.strctProps[i].ProgIDs[j];
//             
//             if (id1 != id2)
//             {
//               m2 = snap2.strctProps[id2].TotMass;
//               
//               ratio = m2/m1;
//               
//               if (ratio > 0.2)
//               {
//                 nmajor++;
//                 tmmajor += m2;
//               }
//               else
//               {
//                 if (ratio > 0.001)
//                 {
//                   nminor++;
//                   tmminor += m2;
//                 }
//                 else
//                 {
//                   nsmooth++;
//                   tmsmooth += m2;
//                 }
//               }
//             }
//             
//             if (ratio > maxratio)
//               maxratio = ratio;
//           }
//           fprintf (fmergers, "%e   %e   %e    %d   %e   %d   %e   %d  %e\n", snap1.strctProps[i].TotMass, m1, maxratio, \
//                    nmajor, tmmajor, nminor, tmminor, nsmooth, tmsmooth);
//         }
//       }
//     }
//     fclose (fmergers);
          
//     int gal1, gal2;
//     double * pos1, * pos2;
//     double rmax1, rmax2;
//     double mass1, mass2;
//     double rsize1, rsize2;
//     double dx;
//     double ratio;
//     double larger;
// 
// 
//     FILE * fpairs  = fopen ("pairs.dat","w");
//     
//     for (i = 1; i <= snap1.nstruct; i++)
//     {
//       if (snap1.strctProps[i].Type == 7)
//       {
//         for (j = 0; j < snap1.strctProps[i].NumSubs; j++)
//         {
//           for (k = j+1; k < snap1.strctProps[i].NumSubs; k++)
//           {
//             gal1 = snap1.strctProps[i].SubIDs[j];
//             gal2 = snap1.strctProps[i].SubIDs[k];
//             
//             mass1 = snap1.strctProps[gal1].TotMass; 
//             mass2 = snap1.strctProps[gal2].TotMass; 
//             
//             pos1 = snap1.strctProps[gal1].Pos;
//             pos2 = snap1.strctProps[gal2].Pos;
//             
//             rmax1 = snap1.strctProps[gal1].Rvmax;
//             rmax2 = snap1.strctProps[gal2].Rvmax;
//             
//             rsize1 = snap1.strctProps[gal1].Rsize;
//             rsize2 = snap1.strctProps[gal2].Rsize;
//             
//             dx = (pos1[0] - pos2[0]) * (pos1[0] - pos2[0]) + \
//                  (pos1[1] - pos2[1]) * (pos1[1] - pos2[1]) + \
//                  (pos1[2] - pos2[2]) * (pos1[2] - pos2[2]);
//             
//             if (mass1 > mass2)
//             {
//               larger = mass1;
//               ratio  = mass2 / mass1;
//             }
//             else
//             {
//               larger = mass2;
//               ratio  = mass1 / mass2;
//             }
//             
//             if (dx < rmax1*rmax1 || dx < rmax2*rmax2)
//               fprintf (fpairs, "%e  %e  %e  %d\n", larger, sqrt(dx), ratio, 0);
//             else
//             {
//               if (dx < rsize1*rsize1 || dx < rsize2*rsize2)
//                 fprintf (fpairs, "%e  %e  %e  %d\n", larger, sqrt(dx), ratio, 1);
//             }
//           }
//         }
//       }
//     }
// 
//     fclose (fpairs);


    
//     free_stfOutput (&snap2);
    
    
      
    
//     load_particles (format, nstruct, strctProps)
//     
//     
//     
    //
    //  Load extended output to extract particles of structures
    //
//     int dummystruct;
//     int dummyindex;
    
    
//     P = malloc ((snap1.nstruct+1) * sizeof(struct pgadget *));
//     for (i = 1; i <= snap1.nstruct; i++)
//     {
//       if (!(P[i] = malloc (snap1.strctProps[i].NumPart * sizeof(struct pgadget))))
//       {
//         printf ("Cannot allocate memory for particle information for structure %d\n", i);
//         printf ("Exiting\n");
//         exit (0);
//       }
//     }
//     
//     int gadget = 0;
//     int ramses = 0;
//     
//     if ( (strcmp(format, "gadget") == 0) )
//       gadget = 1;
//     if ( (strcmp(format, "ramses") == 0) )
//       ramses = 1;
//     
//     int * acum = (int *) malloc ((snap1.nstruct+1) * sizeof(int));
//     for (i = 1; i <= snap1.nstruct; i++)
//       acum[i] = 0;
    
//     int ninextended;
//     if (ramses)
//     {
//       for (i = 0; i < NUMFILES; i++)
//       {
//         ninextended = load_stf_extended_output (stfprefix, i);
//         if (ninextended)
//         {
//           read_ramses_snapshot (directory, galprefix, i);
//           for (j = 0; j < ninextended; j++)
//           {
//             dummystruct = extended_IdStruct[j];
//             dummyindex  = extended_oIndex[j];
//             
//             P[dummystruct][acum[dummystruct]].Pos[0] = ramses_pos[0][dummyindex];
//             P[dummystruct][acum[dummystruct]].Pos[1] = ramses_pos[1][dummyindex];
//             P[dummystruct][acum[dummystruct]].Pos[2] = ramses_pos[2][dummyindex];
//             P[dummystruct][acum[dummystruct]].Vel[0] = ramses_vel[0][dummyindex];
//             P[dummystruct][acum[dummystruct]].Vel[1] = ramses_vel[1][dummyindex];
//             P[dummystruct][acum[dummystruct]].Vel[2] = ramses_vel[2][dummyindex];
//             P[dummystruct][acum[dummystruct]].Mass   = ramses_mass  [dummyindex];
//             P[dummystruct][acum[dummystruct]].Id     = ramses_id    [dummyindex];
//             P[dummystruct][acum[dummystruct]].Age    = ramses_age   [dummyindex];
//                         
//             P[dummystruct][acum[dummystruct]].Type = 4;
//             
//             acum[dummystruct]++;
//           }
//           
//           free_extended_arrays ();
//           free_ramses_arrays ();
//         }
//       }
      
//       double boxlen;
//       double time;
//       double aexp;
//       double H0;
//       double Omega_m;
//       double Omega_l;
//       double Omega_b;
//       char name [NAME_LENGTH];
//       char dummys [NAME_LENGTH];
//       double t;
//       
//       sprintf (name, "%s/info_%s.txt", directory, galprefix);  
//       FILE * ff = fopen(name,"r");
//       for (i = 0; i < 7; i++)
//         fgets(buffer, 100, ff);
//       fgets(buffer, 100, ff);   sscanf (buffer, "%s %s %lf", dummys, dummys, &boxlen);
//       fgets(buffer, 100, ff);   sscanf (buffer, "%s %s %lf", dummys, dummys, &time);
//       fgets(buffer, 100, ff);   sscanf (buffer, "%s %s %lf", dummys, dummys, &aexp);
//       fgets(buffer, 100, ff);   sscanf (buffer, "%s %s %lf", dummys, dummys, &H0);
//       fgets(buffer, 100, ff);   sscanf (buffer, "%s %s %lf", dummys, dummys, &Omega_m);
//       fgets(buffer, 100, ff);   sscanf (buffer, "%s %s %lf", dummys, dummys, &Omega_l);
//       fclose (ff);
//   
//       double * axp_frw;
//       double * hexp_frw;
//       double * tau_frw;
//       double * t_frw;
//       int      n_frw = 1000;
//       double   time_tot = friedman(Omega_m, Omega_l, 0.0, 1e-6, 1e-3, &axp_frw, &hexp_frw, &tau_frw, &t_frw, n_frw);
// 
//       // Find neighbouring conformal time
//       i = 1;
//       while (tau_frw[i] > time  && i < n_frw)
//         i = i+1;
//       
//       // Interpolate time
//       double time_simu = t_frw[i]   * (time - tau_frw[i-1]) / (tau_frw[i]   - tau_frw[i-1]) + \
//                         t_frw[i-1] * (time - tau_frw[i])   / (tau_frw[i-1] - tau_frw[i]);
//                       
//       printf ("Time simu    %lf\n", (time_tot + time_simu) / (H0*1e5/3.08e24) / (365*24*3600*1e9));
//       printf ("Hubble time  %lf\n", time_tot / (H0*1e5/3.08e24) / (365*24*3600*1e9));
//       
//       printf ("i               %d\n", i);
//       printf ("time            %e\n", time);
//       printf ("time_tot        %e\n", time_tot);
//       printf ("time_simu       %e\n", time_simu);
//       printf ("t_tot + t_simu  %e\n", time_tot + time_simu);
//   
//       for (i = 1; i <= snap1.nstruct; i++)
//       {
//         for (j = 0; j < snap1.strctProps[i].NumPart; j++)
//         {
//           k = 1;
//           while (tau_frw[k] > P[i][j].Age  && k < n_frw)
//             k++;
//           
//           t = t_frw[k]   * (P[i][j].Age - tau_frw[k-1]) / (tau_frw[k]   - tau_frw[k-1]) + \
//               t_frw[k-1] * (P[i][j].Age - tau_frw[k])   / (tau_frw[k-1] - tau_frw[k]);
//           P[i][j].Age = (time_simu - t) / (H0*1e5/3.08e24) / (365*24*3600.0);
//         }
//       }
//       free (axp_frw);
//       free (hexp_frw);
//       free (tau_frw);
//       free (t_frw);
//     }
//     char tmpbuffer[100];
    
    
//     if (gadget)
//     {
//       for (i = 0; i < NUMFILES; i++)
//       {
//         if (NUMFILES == 1)
//           sprintf (tmpbuffer, "%s", galprefix);
//         else
//           sprintf (tmpbuffer, "%s.%d", galprefix, i);
//           
//         ninextended = 0;
//         ninextended = load_stf_extended_output (stfprefix, i);
//         
//         if (ninextended)
//         {
//           read_gadget_snapshot (tmpbuffer);
//           for (j = 0; j < ninextended; j++)
//           {
//             dummystruct = extended_IdStruct[j];
//             dummyindex  = extended_oIndex[j];
//             
//             P[dummystruct][acum[dummystruct]].Pos[0] = p[dummyindex].Pos[0];
//             P[dummystruct][acum[dummystruct]].Pos[1] = p[dummyindex].Pos[1];
//             P[dummystruct][acum[dummystruct]].Pos[2] = p[dummyindex].Pos[2];
//                                                       
//             P[dummystruct][acum[dummystruct]].Vel[0] = p[dummyindex].Vel[0];
//             P[dummystruct][acum[dummystruct]].Vel[1] = p[dummyindex].Vel[1];
//             P[dummystruct][acum[dummystruct]].Vel[2] = p[dummyindex].Vel[2];
//             
//             P[dummystruct][acum[dummystruct]].Mass   = p[dummyindex].Mass;
//             P[dummystruct][acum[dummystruct]].Id     = p[dummyindex].Id;
//                         
//             P[dummystruct][acum[dummystruct]].Type = 2;
//             
//             acum[dummystruct]++;
//           }
//           free (p);
//           free (ID);
//           free_extended_arrays ();
//         }
//       }
//     }
    

//     char bfr [NAME_LENGTH];
//     sprintf (bfr, "sSFR_%s.dat", galprefix);
//     FILE * fgal = fopen (bfr, "w");
//     
//     double Mnew;
//     int nage;
//     double Agesum;
//     double mtot;
//     double Dt = 100e6;
//     for (i = 1; i <= snap1.nstruct; i++)
//       if (snap1.strctProps[i].Type != 7)
//       {
//         Mnew = 0;
//         mtot = 0;
//         nage = 0;
//         Agesum = 0;
//         for (j = 0; j < snap1.strctProps[i].NumPart; j++)
//         {
//           if (P[i][j].Age < Dt && P[i][j].Age > 0.0)
// //           if (P[i][j].Age < Dt)
//           {
//             Mnew += P[i][j].Mass;
//             nage++;
//             Agesum += P[i][j].Age;
//           }
//           mtot += P[i][j].Mass;
//         }
//         if (nage)
//           fprintf (fgal, "%d   %e   %e   %e   %e\n", i, snap1.strctProps[i].TotMass, mtot, Mnew/Dt, Agesum/(double)nage);
//       }   
//         
//     fclose (fgal);

    
/*    
    for (i = 1; i <= snap1.nstruct; i++)
      if (snap1.strctProps[i].Type != 7)
      {
        Mnew = 0;
        for (j = 0; j < snap1.strctProps[i].NumPart; j++)
          Mnew += P[i][j].Mass;
        printf ("%d  %e   %e\n", i, snap1.strctProps[i].TotMass, Mnew);
      }   */
    
//     for (i = 1; i <= snap1.nstruct; i++)
//       free (P[i]);
//     free (P);
    
    free_stfOutput (&snap1);
    
  return (0);
}


//---------- Read Propertes File ----------//

void read_properties_file (struct stfOutput * tmpstf)
{
  int    i, j, k;
  
  int    dummyi;
  long   dummyl;
  float  dummyf;
  double dummyd;  
  
  FILE * f;
  char   propts_fname [NAME_LENGTH];
  int    mystructs;

  //
  // Open properties file to read total number of structures
  // and number of processors if stf was run with MPI
  //
  sprintf (propts_fname, "%s.properties", tmpstf->prefix);
  if ((f = fopen (propts_fname, "r")) == NULL)
  {
    sprintf (propts_fname, "%s.properties.0", tmpstf->prefix);
    if ((f = fopen (propts_fname, "r")) == NULL)
    {
      printf ("ERROR: Cannot open file  %s\n", propts_fname);
      exit (0);
    }
  }
  fgets  (longbuffer, NAME_LENGTH, f);  sscanf (longbuffer, "%d  %d", &dummyi, &(tmpstf->nprocs));  
  fgets  (longbuffer, NAME_LENGTH, f);  sscanf (longbuffer, "%d  %d", &dummyi, &(tmpstf->nstruct));  
  fclose (f);

  //
  // Allocate memory for structure properties
  //
  if ( (tmpstf->strctProps = (struct objProps *) malloc ((tmpstf->nstruct+1) * sizeof(struct objProps))) == NULL)
  {
    printf ("Couldn't allocate memory.\n");
    printf ("Exiting...\n");
    exit (0);
  }
  else
    tmpstf->iprops = 1;
  
  for (i = 1; i <= tmpstf->nstruct; i++)
  {
    tmpstf->strctProps[i].SubIDs    = NULL;
    tmpstf->strctProps[i].ProgIDs   = NULL;
    tmpstf->strctProps[i].ProgMrrts = NULL;
  }

  //
  // Read file(s) and store desired structure properties
  //
  int offst = 1;

  for (i = 0; i < tmpstf->nprocs; i++)
  {
    if (tmpstf->nprocs == 1)
      sprintf (propts_fname, "%s.properties", tmpstf->prefix);
    else
      sprintf (propts_fname, "%s.properties.%d", tmpstf->prefix, i);
    
//     printf ("Openning file  %s  \n", propts_fname);
    
    if ((f = fopen (propts_fname, "r")) == NULL)
    {
      printf ("ERROR: Cannot open file  %s", propts_fname);
      exit (0);
    }
    
    fgets  (longbuffer, NAME_LENGTH, f);
    fgets  (longbuffer, NAME_LENGTH, f);
    sscanf (longbuffer, "%d  %d", &mystructs, &dummyi);
    fgets  (longbuffer, 3000, f);
          
    for (j = 0; j < mystructs; j++)
    {
      fgets (longbuffer, 3000, f);
      sscanf (longbuffer, "%d  %d  %d  %d  %d  %d  %d  %lf  %lf  %lf  %lf  %lf  %lf  %lf                \
              %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf            \
              %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf            \
              %lf  %lf",                                                                                \
              &(tmpstf->strctProps[j+offst].ID), &dummyi, &(tmpstf->strctProps[j+offst].DirectHostID),  \
              &(tmpstf->strctProps[j+offst].HostID), &(tmpstf->strctProps[j+offst].NumSubs),            \
              &(tmpstf->strctProps[j+offst].Type), &(tmpstf->strctProps[j+offst].NumPart), &dummyd,     \
              &(tmpstf->strctProps[j+offst].Pos[0]), &(tmpstf->strctProps[j+offst].Pos[1]),             \
              &(tmpstf->strctProps[j+offst].Pos[2]), &dummyd, &dummyd, &dummyd,                         \
              &(tmpstf->strctProps[j+offst].Vel[0]), &(tmpstf->strctProps[j+offst].Vel[1]),             \
              &(tmpstf->strctProps[j+offst].Vel[2]), &dummyd, &dummyd, &dummyd,                         \
              &(tmpstf->strctProps[j+offst].TotMass), &dummyd, &dummyd, &dummyd, &dummyd,               \
              &(tmpstf->strctProps[j+offst].Efrac), &dummyd, &(tmpstf->strctProps[j+offst].Rsize),      \
              &dummyd, &dummyd, &dummyd, &(tmpstf->strctProps[j+offst].RHalfMass),                      \
              &(tmpstf->strctProps[j+offst].Rvmax), &(tmpstf->strctProps[j+offst].Vmax),                \
              &(tmpstf->strctProps[j+offst].Vdisp), &dummyd, &dummyd, &dummyd, &dummyd, &dummyd,        \
              &dummyd, &dummyd, &dummyd, &dummyd, &(tmpstf->strctProps[j+offst].Lambda),                \
              &(tmpstf->strctProps[j+offst].L[0]), &(tmpstf->strctProps[j+offst].L[1]),                 \
              &(tmpstf->strctProps[j+offst].L[2])                                                       \
            );
    }    
    offst += mystructs;
    fclose(f);    
  } 
}



//---------- Read Gadget Snapshot ---------//
void read_gadget_snapshot(char * snapshot)
{
  FILE *fd;
  int k, dummy, ntot_withmasses;
  int n, pc, pc_new;
#define SKKIP fread(&dummy, sizeof(dummy), 1, fd);

  pc = 0;

  if((fd = fopen(snapshot, "r")) == NULL)
  {
    printf("Couldn't open file\n");
    exit(0);
  }
  
  fflush(stdout);
  
  fread(&dummy, sizeof(dummy), 1, fd);
  fread(&header1, sizeof(header1), 1, fd);
  fread(&dummy, sizeof(dummy), 1, fd);  
  
  for(k = 0, NumPart = 0, ntot_withmasses = 0; k < 6; k++)
    NumPart += header1.npart[k];

  for(k = 0, ntot_withmasses = 0; k < 6; k++)
    if(header1.mass[k] == 0)
      ntot_withmasses += header1.npart[k];
  
  printf("Allocating memory...");

  if(!(p = malloc(NumPart * sizeof(struct pgadget))))
    {
      fprintf(stderr, "failed to allocate memory.\n");
      exit(0);
    }
    
  if(!(ID = malloc(NumPart * sizeof(int))))
    {
      fprintf(stderr, "failed to allocate memory.\n");
      exit(0);
    }

  printf("done\n");
  
  SKKIP;

  for(k = 0, pc_new = pc; k < 6; k++)
    for(n = 0; n < header1.npart[k]; n++)
    {
      fread(&p[pc_new].Pos[0], sizeof(float), 3, fd);
      pc_new++;
    }
  
  SKKIP;
  SKKIP;

  for(k = 0, pc_new = pc; k < 6; k++)
    for(n = 0; n < header1.npart[k]; n++)
    {
      fread(&p[pc_new].Vel[0], sizeof(float), 3, fd);
      pc_new++;
    }
  
  SKKIP;
  SKKIP;

  for(k = 0, pc_new = pc; k < 6; k++)
    for(n = 0; n < header1.npart[k]; n++)
    {
      fread(&p[pc_new].Id, sizeof(int), 1, fd);
      ID[pc_new] = p[pc_new].Id;
      pc_new++;
    }
  SKKIP;
    
  if(ntot_withmasses>0)
    SKKIP;
  for(k = 0, pc_new = pc; k < 6; k++)
  {
    for(n=0; n < header1.npart[k]; n++)
    {
      p[pc_new].Type = k;
      if(header1.mass[k] == 0)
	fread(&p[pc_new].Mass, sizeof(float), 1, fd);
      else
	p[pc_new].Mass= header1.mass[k];
      pc_new++;
    }
  }
  if(ntot_withmasses>0)
    SKKIP;

  
  int check;
  for(k = 0, pc_new = pc; k < 6; k++)
  {
    check = 0;
    for(n = 0; n < header1.npart[k]-1; n++)
    {
      if(p[pc_new + 1].Mass == p[pc_new].Mass)
	check++;
      pc_new++;
    }
    if(header1.npart[k] > 0) pc_new++;
    if((check == header1.npart[k]-1) && (header1.mass[k] == 0))
      header1.mass[k] = p[pc_new-1].Mass;
  }
//   printf("NumPart =   %i\n", NumPart);
//   printf("PCNEW =     %i\n", pc_new);
  printf("Successful snapshot reading\n");
  
  fclose(fd);
}

//---------- Write Snapshot File ----------//

void write_snapshot(struct pgadget * pg, int particles, struct io_header_1 header, char * output_name)
{
  FILE * snap_file;
  
  int dummy;
  int k, pc_new, n;
  int pc = 0;
  int ntot_withmasses;
  
  double ref_mass[6];
  int id_ref[6];
  int ** ids;
  
  int i, bob, offset;

  if((snap_file = fopen(output_name,"w")) == NULL)
  {
    printf("Couldn't open file\n");
    exit(0);
  }
    
  //! Initialize
  for(k = 0; k < 6; k++)
  {
    header.npart[k]      = 0;
    header.npartTotal[k] = 0;
    header.mass[k]       = 0;
    ref_mass[k]          = 0;
    id_ref[k]            = 0;
  }
  header.time           = 0.9823335728;
  header.redshift       = 0.0179841427;
  header.flag_sfr       = 0;
  header.flag_feedback  = 0;
  header.flag_cooling   = 0;
  header.num_files      = 1;
  header.BoxSize        = 100000;
  header.Omega0         = 0.2720000148;
  header.OmegaLambda    = 0.7279999852;
  header.HubbleParam    = 0.7040000153;
  
  
  //! Count Particles by type
  for(k = 0; k < particles; k++)
  {
    header.npart[pg[k].Type]++;
    header.npartTotal[pg[k].Type]++;
  }
  
  ids = (int **)malloc(6 * sizeof(int *));
  for(k = 0; k < 6; k++)
   ids[k] = (int *)malloc((header.npart[k]) * sizeof(int));

  //! CHECK MASSES  
  for(n = 0; n < particles; n++)
  {
    header.mass[pg[n].Type] += pg[n].Mass;
    ref_mass[pg[n].Type] = pg[n].Mass;
    
    bob = pg[n].Type;
    
    ids[bob][id_ref[bob]] = n;
    id_ref[bob]++;
  }
  
  for(k = 0; k < 6; k++)
  {
    if((header.mass[k] - ref_mass[k] * header.npart[k] == 0) && header.npart[k] > 0)
      header.mass[k] = ref_mass[k];
    else
      header.mass[k] = 0;    
//     printf("%i\n", k);
//     printf("      Npart  %i\n", header.npart[k]);
//     printf("      Mass   %g\n", header.mass[k]);
  }

//   offset = 0;
//   for(k = 0; k < 6; k++)
//   {
//     if(header.npart[k] > 0)
//       for(i = 0; i < header.npart[k]; i++)
// 	ID[offset + i] = ids[k][i];
//       
//     offset += header.npart[k];
//   }
      
  //! START WRITING  
  fflush(stdout);
  header.num_files = 1;

  dummy = 256;
  fwrite(&dummy, sizeof(dummy), 1, snap_file);
  fwrite(&header, sizeof(header), 1, snap_file);
  fwrite(&dummy, sizeof(dummy), 1, snap_file);

//   for(k = 0; k < 6; k++)
//     printf("\tType %i Particles \t%i\n", k, header.npart[k]);
//   for(k = 0; k < 6; k++)
//     printf("\tType %i Mass      \t%g\n", k, header.mass[k]);
//   printf("\tTime                \t%g\n", header.time);
//   printf("\tRedshift            \t%g\n", header.redshift);
//   printf("\tSFR                 \t%i\n", header.flag_sfr);
//   printf("\tFeedback            \t%i\n", header.flag_feedback);
//   for(k = 0; k < 6; k++)
//     printf("\tType %i npart Total \t%i\n", k, header.npartTotal[k]);
//   printf("\tCooling             \t%i\n", header.flag_cooling);
//   printf("\tnum files           \t%i\n", header.num_files);
//   printf("\tBox Size            \t%g\n", header.BoxSize);
//   printf("\tOmega0              \t%g\n", header.Omega0);
//   printf("\tOmegaL              \t%g\n", header.OmegaLambda);
//   printf("\tHubbleParam         \t%g\n", header.HubbleParam);
//   printf("\tFill bytes          \t%s\n", header.fill);
// 

  //!------ Pos
  dummy = 3 * 4 * particles;
  fwrite(&dummy, sizeof(dummy), 1, snap_file);
  for(k = 0; k < particles; k++)
    fwrite(&pg[k].Pos[0], sizeof(float), 3, snap_file);
  fwrite(&dummy, sizeof(dummy), 1, snap_file);
  
  //!------ Vel
  fwrite(&dummy, sizeof(dummy), 1, snap_file);
  for(k = 0; k < particles; k++)
    fwrite(&pg[k].Vel[0], sizeof(float), 3, snap_file);  
  fwrite(&dummy, sizeof(dummy), 1, snap_file);
  
  //!------- ID
  dummy = 4 * particles;
  fwrite(&dummy, sizeof(dummy), 1, snap_file);
  for(k = 0; k < particles; k++)
    fwrite(&pg[k].Id, sizeof(int), 1, snap_file);
  fwrite(&dummy, sizeof(dummy), 1, snap_file);

   //!------- Mass  
  dummy = 0;
  for(k = 0; k < 6; k++)
    if((header.mass[k] == 0) && (header.npart[k] > 0))
      dummy += header.npart[k];
  dummy *= sizeof(float);
  
  offset = 0;
  if(dummy != 0)  
  {
//     printf("writing mass block\n");
    fwrite(&dummy, sizeof(dummy), 1, snap_file);
    for(k = 0; k < 6; k++)
    {
      if((header.npart[k] > 0) && (header.mass[k] == 0))
      {
	for(i = 0; i < header.npart[k]; i++)
	  fwrite(&pg[k].Mass, sizeof(float), 1, snap_file);
      }
    }
    fwrite(&dummy, sizeof(dummy), 1, snap_file);
  }

//   printf("done writing\n");
 
  //! free memory
  for(k = 0; k < 6; k++)
    if(header.npart[k] > 0)
      free(ids[k]);
  free(ids);
  
  //! Close file
  fclose(snap_file);  
}

// ---  Read Galfile --- //

void read_gal_file (char * filename)
{  
  int     dummy;
  int     i, j;
  FILE *  f;
  
  if ((f = fopen (filename, "r")) == NULL)
  {
    printf ("Can't open file named   %s \n", filename);
    exit(0);
  }

#define SKIP fread (&dummy, sizeof(int), 1, f);  
#define PSKIP printf ("Block-size  %d \n", dummy);
     
  SKIP  fread (&gal_number, sizeof(int),    1, f);  SKIP
  SKIP  fread (&gal_level,  sizeof(int),    1, f);  SKIP  
  SKIP  fread (&gal_mass,   sizeof(double), 1, f);  SKIP
  SKIP  fread (&gal_pos,    sizeof(double), 3, f);  SKIP
  SKIP  fread (&gal_vel,    sizeof(double), 3, f);  SKIP
  SKIP  fread (&gal_ang,    sizeof(double), 3, f);  SKIP
  SKIP  fread (&nlist,      sizeof(int),    1, f);  SKIP
  
  printf ("\n");
  printf ("my_number   %d \n", gal_number);
  printf ("Level       %d \n", gal_level);
  printf ("Mass        %8.5f \n", gal_mass);
  printf ("Pos         %8.5f \t %8.5f \t %8.5f \n", gal_pos[0], gal_pos[1], gal_pos[2]);
  printf ("Vel         %8.5f \t %8.5f \t %8.5f \n", gal_vel[0], gal_vel[1], gal_vel[2]);
  printf ("Ang Mom     %8.5f \t %8.5f \t %8.5f \n", gal_ang[0], gal_ang[1], gal_ang[2]);
  printf ("Nlist       %d \n", nlist);
  printf ("\n");

  if (!(p = malloc (nlist * sizeof(struct pdata))))
  {
    printf ("Cannot allocate memory for particle information\n");
    printf ("Exiting\n");
    exit (0);
  }

  for (j = 0; j < 3; j++)
  {
    SKIP
    if (dummy == 8 * nlist)
      for (i = 0; i < nlist; i++)
      {
	fread (&p[i].Pos[j], sizeof(double), 1, f);
	p[i].Pos[j] = p[i].Pos[j] * 1000; 
      }
    SKIP
  }
  
  for (j = 0; j < 3; j++)
  {
    SKIP
    if (dummy == 8 * nlist)
      for (i = 0; i < nlist; i++)
	fread (&p[i].Vel[j], sizeof(double), 1, f);
    SKIP
  }        

  SKIP
  if (dummy == 8 * nlist)
    for (i = 0; i < nlist; i++)
    {
      fread (&p[i].Mass, sizeof(double), 1, f);
      p[i].Mass *= 10.0;                          // converts to M/10**10 Msun
    }
  SKIP
  
  SKIP
  if (dummy == 4 * nlist)
    for (i = 0; i < nlist; i++)
      fread (&p[i].Id, sizeof(int), 1, f);
  SKIP
  
  SKIP
  if (dummy == 8 * nlist)
    for (i = 0; i < nlist; i++)
      fread (&p[i].Age, sizeof(double), 1, f);
  SKIP

  SKIP
  if (dummy == 8 * nlist)
    for (i = 0; i < nlist; i++)
      fread (&p[i].Metal, sizeof(double), 1, f);
  SKIP
  
  SKIP
  if (dummy == 8 * nlist)
    for (i = 0; i < nlist; i++)
      fread (&p[i].Chem, sizeof(double), 1, f);
  SKIP
  
  SKIP PSKIP   
  fclose(f); 
}


void read_ramses_snapshot(char * dir, char * snapshot, int numfile)
{
  char name[100];
  char buffer[100];
  char dummys[100];
  int i,j;
  
  double boxlen;
  double time;
  double aexp;
  double H0;
  double Omega_m;
  double Omega_l;
  double Omega_b;
  double unit_l;
  double unit_d;
  double unit_v;
  double unit_t;
  double unit_m;
  
  sprintf (name, "%s/info_%s.txt", dir, snapshot);  
  FILE * ff = fopen(name,"r");
  for (i = 0; i < 7; i++)
    fgets(buffer, 100, ff);
  
  fgets(buffer, 100, ff);   sscanf (buffer, "%s %s %lf", dummys, dummys, &boxlen);
  fgets(buffer, 100, ff);   sscanf (buffer, "%s %s %lf", dummys, dummys, &time);
  fgets(buffer, 100, ff);   sscanf (buffer, "%s %s %lf", dummys, dummys, &aexp);
  fgets(buffer, 100, ff);   sscanf (buffer, "%s %s %lf", dummys, dummys, &H0);
  fgets(buffer, 100, ff);   sscanf (buffer, "%s %s %lf", dummys, dummys, &Omega_m);
  fgets(buffer, 100, ff);   sscanf (buffer, "%s %s %lf", dummys, dummys, &Omega_l);
  fgets(buffer, 100, ff);   sscanf (buffer, "%s %s %s ", dummys, dummys, dummys);
  fgets(buffer, 100, ff);   sscanf (buffer, "%s %s %lf", dummys, dummys, &Omega_b);
  fgets(buffer, 100, ff);   sscanf (buffer, "%s %s %lf", dummys, dummys, &unit_l);
  fgets(buffer, 100, ff);   sscanf (buffer, "%s %s %lf", dummys, dummys, &unit_d);
  fgets(buffer, 100, ff);   sscanf (buffer, "%s %s %lf", dummys, dummys, &unit_t);  
  fclose (ff);
  
  unit_m = unit_d * unit_l * unit_l * unit_l;  // in grams
  unit_m = unit_m / 1.989e+33;                 // in solar masses
  
  unit_v = unit_l / unit_t;                    // in cm / s
  unit_v = unit_v / 100000.0;                  // in km / s
  
  unit_l = unit_l / 3.08e+21;
//   printf("DM mass  %g\n", 7.755314938e-10 * unit_m);
  
  sprintf (name, "%s/part_%s.out%05d", dir, snapshot, numfile+1);
//   printf ("Reading particle file   %s\n", name);

  FILE * f = fopen(name,"r");
  
//   fread (ptr, num, nmemb, stream)
  int dummy;
  int info[100];
  
#define SSKIP   dummy=0; fread (&dummy, sizeof(int), 1, f);
#define PSSKIP  printf ("%d  ", dummy);
#define PSSKIP2 printf ("%d\n", dummy);

  int lala;
 
  int size = 0;
  
  int     ncpu;
  int     npart;
  int     seed[4];
  int     nstarTot;
  double  mstarTot;
  double  mstarLst;
  int     nsink;
    
  //!--- Header
  SSKIP  fread(&ncpu,     sizeof(int),    1, f);  SSKIP     
  SSKIP  fread(&ndim,     sizeof(int),    1, f);  SSKIP     
  SSKIP  fread(&npart,    sizeof(int),    1, f);  SSKIP   
  SSKIP  fread(&seed[0],  sizeof(int),    4, f);  SSKIP   
  SSKIP  fread(&nstarTot, sizeof(int),    1, f);  SSKIP   
  SSKIP  fread(&mstarTot, sizeof(double), 1, f);  SSKIP   
  SSKIP  fread(&mstarLst, sizeof(double), 1, f);  SSKIP   
  SSKIP  fread(&nsink,    sizeof(int),    1, f);  SSKIP   

//   printf ("NumProcs        %d\n", ncpu);
//   printf ("Num Dims        %d\n", ndim);
//   printf ("Npart           %d\n", npart);
//   for (i = 0; i < 4; i++)
//     printf ("LocalSeed[%d]    %d\n", i, seed[i]);    
//   printf ("NstarTot        %d\n", nstarTot);
//   printf ("Mstar_tot       %g\n", mstarTot);
//   printf ("Mstar_lost      %g\n", mstarLst);
//   printf ("NumSink         %d\n", nsink);
  
  //!--- Allocate for DM and Stars  
  ramses_pos = (double **) malloc (ndim * sizeof(double *));
  ramses_vel = (double **) malloc (ndim * sizeof(double *));
  ramses_met = (double **) malloc (11   * sizeof(double *));
  
  for (i = 0; i < ndim; i++) ramses_pos[i] = (double *) malloc (npart * sizeof(double));
  for (i = 0; i < ndim; i++) ramses_vel[i] = (double *) malloc (npart * sizeof(double));
  for (i = 0; i < 11;   i++) ramses_met[i] = (double *) malloc (npart * sizeof(double));
  
  ramses_age  = (double *) malloc (npart * sizeof(double));
  ramses_mass = (double *) malloc (npart * sizeof(double));
  ramses_id   = (int    *) malloc (npart * sizeof(int));
  ramses_lvl  = (int    *) malloc (npart * sizeof(int));

  //--- Pos
  for (i = 0; i < ndim; i++)
  {
    SSKIP  fread(&ramses_pos[i][0], sizeof(double), npart, f);  SSKIP
  }
  //--- Vel
  for (i = 0; i < ndim; i++) 
  {
    SSKIP  fread(&ramses_vel[i][0], sizeof(double), npart, f);  SSKIP
  }
  //--- Mass
  SSKIP  fread(&ramses_mass[0], sizeof(double), npart, f);  SSKIP
  
  //--- Id
  SSKIP  fread(&ramses_id[0], sizeof(int), npart, f);  SSKIP

  //--- Level
  SKIP  fread(&ramses_lvl[0], sizeof(int), npart, f);  SKIP
  
  //--- Birth Epoch
  SKIP  fread(&ramses_age[0], sizeof(double), npart, f);  SKIP
  
  //--- Metallicity if ((STAR || SINK) && (METAL))
  for (i = 0; i < 11; i++)
  {
    SKIP  fread(&ramses_met[i][0], sizeof(double), npart, f);  SKIP
  }
  fclose (f);

  
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
//                      t_frw[i-1] * (time - tau_frw[i])   / (tau_frw[i-1] - tau_frw[i]);
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
  // Convert to human readable units
  //
  double t;
  for (i = 0; i < npart; i++)
  {
    ramses_pos[0][i] *= unit_l;
    ramses_pos[1][i] *= unit_l;
    ramses_pos[2][i] *= unit_l;
    
    ramses_vel[0][i] *= unit_v;
    ramses_vel[1][i] *= unit_v;
    ramses_vel[2][i] *= unit_v;
    
    ramses_mass[i]   *= unit_m;
    
//     if (ramses_age[i] != 0)
//     {    
//       j = 1;
//       while (tau_frw[j] > ramses_age[i]  && j < n_frw)
//         j++;
//       
//       printf ("%e  ", ramses_age[i]);
//       
//       t = t_frw[j]   * (ramses_age[i] - tau_frw[j-1]) / (tau_frw[j]   - tau_frw[j-1]) + \
//           t_frw[j-1] * (ramses_age[i] - tau_frw[j])   / (tau_frw[j-1] - tau_frw[j]);
//       ramses_age[i] = (time_simu - t) / (H0*1e5/3.08e24) / (365*24*3600.0);
//       printf ("%e  %e\n", t, ramses_age[i]);
//     }
  }
  
//   free (axp_frw);
//   free (hexp_frw);
//   free (tau_frw);
//   free (t_frw);
  
  
}


int read_stf_filesofgroup (char * prefix, int strct_id, int ** files_of_strct)
{
  int i, j, k;
  FILE * f;
  
  char buffer [NAME_LENGTH];
  
  sprintf (buffer, "%s.filesofgroup", prefix);
  
  int tmpid;
  int nfiles;
  
  f = fopen (buffer, "r");
  do
  {
    fgets (buffer, NAME_LENGTH, f);
    sscanf (buffer, "%d  %d", &tmpid, &nfiles);
    fgets (buffer, NAME_LENGTH, f);
    get_n_num_from_string (buffer, nfiles, files_of_strct);
  }
  while (tmpid != strct_id);
  
  fclose (f);
  
  return nfiles;
}


int get_n_num_from_string (char * strng, int n_num, int ** nums)
{
  int  i, j, k;
  int  ndigits = 10;
  char tmpbuff [ndigits];
  int  length = strlen (strng);
  
  for (i = 0; i < ndigits; i++)
    tmpbuff[i] = '\0';
  
  for (i = 0, k = 0; i < length; i++)
    if (strng[i] == ' ')
      k++;
    
  if (k != n_num)
  {
    printf ("ERROR: there should be %d numbers, but there are only %d\n", n_num, k);
    printf ("Exiting...\n");
  }
  
  *nums = (int *) malloc (k * sizeof(int));
  
  for (i = 0, j = 0, k = 0; i < length; i++)
  {
    if (strng[i] != ' ')
    {
      tmpbuff[j] = strng[i];
      tmpbuff[j+1] = '\0';
      j++;
    }
    else
    {
      nums[0][k] = atoi (tmpbuff);
      k++;
      j = 0;
    }
  }

  
  return 0;
  
}


int load_treefrog (char * tffile, int strct_id, int ** prog_ids, float ** prog_mrrts)
{
  int i, j, k;
  FILE * f;
  
  int tmpid;
  int nprogs;
  
  char buffer [LONG_LENGTH];
  
  f = fopen (tffile, "r");
  
  fgets (buffer, LONG_LENGTH, f);
  fgets (buffer, LONG_LENGTH, f);
  fgets (buffer, LONG_LENGTH, f);
  fgets (buffer, LONG_LENGTH, f);
  
  do
  {
    fgets (buffer, LONG_LENGTH, f);
    sscanf (buffer, "%d  %d", &tmpid, &nprogs);
    
    if (tmpid != strct_id)
      for (i = 0; i < nprogs; i++)
        fgets (buffer, LONG_LENGTH, f);
  }
  while (tmpid != strct_id);
  
  
  *prog_ids   = (int *)   malloc (nprogs * sizeof(int));
  *prog_mrrts = (float *) malloc (nprogs * sizeof(float));
  
  for (i = 0; i < nprogs; i++)
  {
    fgets (buffer, LONG_LENGTH, f);
    sscanf (buffer, "%d  %f", &prog_ids[0][i], &prog_mrrts[0][i]);    
  }  
  fclose (f);
  
  return nprogs;
}


int load_stf_extended_output (char * prefix, int filenum)
{
    FILE * f;
    char buffer[NAME_LENGTH];
    sprintf(buffer, "%s.extended.%d", prefix, filenum);
    int nparts = 0;
    int i;
    
    if ((f = fopen(buffer, "r")) == NULL)
      return 0;
    
    while (fgets(buffer, NAME_LENGTH, f) != NULL)
      nparts++;
    rewind(f);
    
//     printf("nparts %d\n", nparts);
    
    extended_oIndex   = (int *) malloc (nparts * sizeof(int));
    extended_IdStruct = (int *) malloc (nparts * sizeof(int));
    extended_IdHost   = (int *) malloc (nparts * sizeof(int));
    extended_IdIGM    = (int *) malloc (nparts * sizeof(int));
    
    for (i = 0; i < nparts; i++)
    {
      fgets(buffer, NAME_LENGTH, f);
      sscanf(buffer, "%d  %d  %d  %d  ", &extended_oIndex[i], &extended_IdStruct[i], &extended_IdHost[i], &extended_IdIGM[i]);
    }
    
    fclose (f);
    
    return nparts;
}

void free_extended_arrays (void)
{
  free (extended_oIndex);
  free (extended_IdStruct);
  free (extended_IdHost);
  free (extended_IdIGM);
}

void free_ramses_arrays (void)
{
  int i;
  
  for (i = 0; i < ndim; i++)
  {
    free (ramses_pos[i]);
    free (ramses_vel[i]);
  }
  for (i = 0; i < 11; i++)
    free (ramses_met[i]);
  
  free (ramses_pos);
  free (ramses_vel);
  free (ramses_mass);
  free (ramses_id); 
  free (ramses_met);
  free (ramses_age);
  free (ramses_lvl);
}


void init_stfOutput (struct stfOutput * tmpstf)
{
  tmpstf->nprocs     = 0;
  tmpstf->nstruct    = 0;
  tmpstf->iprops     = 0;
  tmpstf->iparts     = 0;
  tmpstf->strctProps = NULL;
  tmpstf->strctParts = NULL;
}
    
void free_stfOutput (struct stfOutput * tmpstf)
{
  int i;
  
  if (tmpstf->iprops)
  {
    for (i = 1; i <= tmpstf->nstruct; i++)
    {
      if (tmpstf->strctProps[i].SubIDs != NULL)
        free (tmpstf->strctProps[i].SubIDs);
      
      if (tmpstf->strctProps[i].ProgIDs != NULL)
        free (tmpstf->strctProps[i].ProgIDs);
      
      if (tmpstf->strctProps[i].ProgMrrts != NULL)
        free (tmpstf->strctProps[i].ProgMrrts);
    }
  }
  
  if (tmpstf->iparts)
  {
    for (i = 1; i <= tmpstf->nstruct; i++)
      free (tmpstf->strctParts[i]);
    free (tmpstf->strctParts);
  }
}
    

    
void fill_SubIDS (struct stfOutput * tmpstf)
{
  int i, j;
  int bob;
  int tmp;

  
  for (i = 1; i <= tmpstf->nstruct; i++)
    if (tmpstf->strctProps[i].HostID == -1) 
    {
      bob = tmpstf->strctProps[i].NumSubs;
      tmpstf->strctProps[i].SubIDs = (int *) malloc (bob * sizeof(int));
      tmpstf->strctProps[i].dummy = 0;
      for (j = 0; j < bob; j++)
        tmpstf->strctProps[i].SubIDs[j] = 0;
    }
    
  for (i = 1; i <= tmpstf->nstruct; i++)
  {
    if (tmpstf->strctProps[i].HostID != -1) 
    {
      bob = i;
      while (tmpstf->strctProps[bob].HostID != -1)
      {
        bob = tmpstf->strctProps[bob].DirectHostID;
      }
      tmp = tmpstf->strctProps[bob].dummy++; 
      tmpstf->strctProps[bob].SubIDs[tmp] = i;
    }
  } 
}


void fill_ProgIDs (struct stfOutput * tmpstf, char * tffile)
{
  int i, j;
  int bob;
  int tmpid ;
  int nprogs; 
  
  FILE * f = fopen (tffile, "r");
  
  fgets (buffer, LONG_LENGTH, f);
  fgets (buffer, LONG_LENGTH, f);
  fgets (buffer, LONG_LENGTH, f);
  fgets (buffer, LONG_LENGTH, f);
  
  for (i = 1; i <= tmpstf->nstruct; i++)
  {
    fgets (buffer, LONG_LENGTH, f);
    sscanf (buffer, "%d  %d", &tmpid, &nprogs);
    
    tmpstf->strctProps[i].NumProg = nprogs;
    if (nprogs > 0)
    {
      tmpstf->strctProps[i].ProgIDs   = (int *)    malloc (nprogs * sizeof(int));
      tmpstf->strctProps[i].ProgMrrts = (double *) malloc (nprogs * sizeof(double));
      for (j = 0; j < nprogs; j++)
      {
        fgets  (buffer, LONG_LENGTH, f);
        sscanf (buffer, "%d  %lf", &(tmpstf->strctProps[i].ProgIDs[j]), &(tmpstf->strctProps[i].ProgMrrts[j]));
      }
    }
  }
  
  fclose (f);
}



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


double dadtau (double axp_tau, double Omega0, double OmegaL, double OmegaK)
{
  return sqrt(axp_tau*axp_tau*axp_tau * (Omega0 + OmegaL*axp_tau*axp_tau*axp_tau + OmegaK*axp_tau));
}

double dadt (double axp_t, double Omega0, double OmegaL, double OmegaK)
{
  return sqrt((1.0/axp_t) * (Omega0 + OmegaL*axp_t*axp_t*axp_t + OmegaK*axp_t));
}



//     FILE * ff;
//     int PID, STID, HSTID, IGMID;
//     
//     printf("Distributing particles\n");
//     for (i = 0; i < NUMFILES; i++)
//     {
//       sprintf (buffer, "%s.%d", galprefix, i);
//       read_gadget_snapshot (buffer);
//       sprintf (buffer, "%s.extended.%d", stfprefix, i);
//       ff = fopen (buffer, "r");
//       for (j = 0; j < NumPart; j++)
//       {
// 	fgets  (longbuffer, NAME_LENGTH, ff);
// 	sscanf (longbuffer, "%d %d %d %d", &PID, &STID, &HSTID, &IGMID);
// 	if (STID > 0)
// 	{
// 	  P[STID-1][acum[STID-1]].Id = p[j].Id;
// 	  P[STID-1][acum[STID-1]].Pos[0] = p[j].Pos[0];
// 	  P[STID-1][acum[STID-1]].Pos[1] = p[j].Pos[1];
// 	  P[STID-1][acum[STID-1]].Pos[2] = p[j].Pos[2];
// 	  P[STID-1][acum[STID-1]].Vel[0] = p[j].Vel[0];
// 	  P[STID-1][acum[STID-1]].Vel[1] = p[j].Vel[1];
// 	  P[STID-1][acum[STID-1]].Vel[2] = p[j].Vel[2];
// 	  P[STID-1][acum[STID-1]].Mass = p[j].Mass;
// 	  P[STID-1][acum[STID-1]].Type = p[j].Type;
// 	  acum[STID-1]++;
// 	}
//       }  
//       fclose(ff);
//       free (p);
//     }
//     printf("Particles Distributed\n");
    
//     //
//     // Write Gadget-formatted files
//     //
//     for (i = 0; i < nstruct; i++)
//     {
//       printf ("struct %d nparts %d acum %d \n", i, nparts[i], acum[i]);
//       if (outType == 0)
//         sprintf (output_fname, "%s_%03d.gdt_%03d", outprefix, hostid[i], idins[i]-1);
//       if (outType == 1)
//         sprintf (output_fname, "%s_%03d", outprefix, i+1);
//       write_snapshot(P[i], acum[i], header1, output_fname);
//     }
