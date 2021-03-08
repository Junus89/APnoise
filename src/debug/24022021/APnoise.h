#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<complex.h>
#include<time.h>
#include<string.h>

#define PI      3.14159265358979323846
#define TWOPI   (2.0*PI)
#define DEBUG 1

/* User defined data types */
  typedef struct{
    double x[2];
  }d2_t;

  typedef struct{
   double minVal;
   int    minValIdx;
  }minValLoc;

  typedef struct{
    double complex x[2];
  }dc2_t;

  typedef struct{
    double x[3];
  }d3_t;

  typedef struct{
    double complex x[3];
  }dc3_t;

  typedef struct{
    double x[4];
  }d4_t;

  typedef struct{
    double x[5];
  }d5_t;

  typedef struct{
    double x[4];
    char   FoilFN[128];
  }Airfoil_t;


  typedef struct{
    double complex x[5];
  }dc5_t;

  typedef struct{
    double x[6];
  }d6_t;

  typedef struct{
    double complex x[6];
  }dc6_t;


typedef struct{
  char caseName[256];
  int runID;
}caseData;
typedef struct{
  double  D;
  int     NB;
  int     RPM;
  double  Mach;
  double  CT;
  double  eta;
  double  omega;
  double  Thrust;
  double  Torque;
  double  Re;
  double  bladeAngle;
  double  bladeChord;
}PropellerData;

typedef struct{
  double  alt;
  double  T;
  double  p;
  double  c0;
  double  rho;
}AtmData;

typedef struct{
  int     thetaNum;
  int     phiNum;
  int     NFNum;
  int     FFNum;
  double  *theta;
  double  *phi;
  d3_t    *NFcoords;/* near field */
  d3_t    *FFcoords;/* far field */
}ObserverData;
  
enum runtype{
  GutinDeming=0,\
  GarrickWatkins=1,\
  BarryMagliozzi=2,\
  Hanson=3 
};





/* Utility functions */

  void APnoise_error(char error_text[],\
                     const char* routine);
  void APnoise_warning(char error_text[],\
                      const char* routine);

  /*------------------------------------------------------------------*/
  /* Functions/Routines for memory allocation for 1D-4D arrays*/

  /*--------------------------------1D--------------------------------*/
  double *dvector(long nl,long nh);

  void free_dvector(double *inputAr,\
                    long nl,long nh);

  /* complex dvector */
  double complex *dcvector(long nl,long nh);

  void free_dcvector(double complex *inputAr,\
                     long nl,long nh);

  /*--------------------------------2D--------------------------------*/

  double **dmatrix(long nrl,long nrh, \
                   long ncl,long nch);
  /* */
  void free_dmatrix(double **inputAr,\
                    long nrl,long nrh,\
                    long ncl,long nch);

  d2_t **d2_tmatrix(long nrl,\
                    long nrh,\
                    long ncl,\
                    long nch);
  /* */
  void free_d2_tmatrix(d2_t **inputAr,\
                       long nrl,\
                       long nrh,\
                       long ncl,\
                       long nch);
/*------------------------------3D------------------------------------*/

  double ***d3dtensor(long nrl,long nrh,\
                      long ncl,long nch,\
                      long ndl,long ndh);

  void free_d3dtensor(double ***inputAr,\
                      long nrl,long nrh,\
                      long ncl,long nch,\
                      long ndl,long ndh);

  /* */
  /* complex dmatrix */
  double complex **dcmatrix(long nrl,long nrh,\
                            long ncl,long nch);

  void free_dcmatrix(double complex **inputAr,\
                     long nrl,long nrh,\
                     long ncl,long nch);


  /* Mem functions for user defined data types */
  d2_t *d2_tvector(long nl,long nh);
  void free_d2_tvector(d2_t *inputAr,long nl,long nh);

  d3_t *d3_tvector(long nl,long nh);
  void free_d3_tvector(d3_t *inputAr,long nl,long nh);
/*--------------------------- File reader utility functions ----------*/
  d3_t *freader3col(char* fn,int* n);

  d4_t *polarDataReader(char* fn,int* n);

  Airfoil_t *AirfoilDBreader(char* fn,int* n);


/* --------------------------- Math utility functions ----------------*/
  d2_t blii1d(d2_t x,\
              d2_t y,
              double xi);

  d3_t bcui1d(d4_t x,\
              d4_t f,\
              double xi);

  d2_t blii1d(d2_t x,\
              d2_t y,
              double xi);

  int brcket(int n,\
             double x[],\
             double xi);

  d3_t interp1D(int nx,\
                double tabx[],\
                double tabc[],\
                double xi);

  double *diff(double *vec,\
               int N);

  double *linspace(double start,\
                   double end,\
                   int N);

  minValLoc findMinValLoc(double *arr,\
                          int numArr);

  char* cutoffstr(const char* str, int from , int to);

  int strFind(Airfoil_t *profile,\
              int  N,\
              int  idx,\
              char pattern[128]);

  double integTrapz(double *x, int N);

  double factorial(int n);

  double besselj(int n,double x);

  double get_cpu_time(clock_t start,\
                      clock_t end);


