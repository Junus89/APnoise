#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<complex.h>
#include<time.h>
#include<string.h>

#define PI      3.14159265358979323846
#define TWOPI   (2.0*PI)
#define DEBUG 1

#define MIN(a,b) (((a)<(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))

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
  char Fn[256];
  int NumStations;
  double *r_R;
  double *r;
  double *c_R;
  double *c_D;
  double *c;
  double *twist;
  double *thickns;
  double thicknsRatio;
  double *lambdaManufac;
  double R;
  double NB;
  double hub_r;
  double RPM;
  double U;
}BladeGeomData;

typedef struct{
  char Fn[256];
  int NumStations;
  double *Cl;
  double *Cd;
  double *AoA;
  double *Veff;
  double *Ct;
  double *Cq;
}BladeAeroData;

typedef struct{
  int    NumStations;        /* number of stations */
  double *tb;                /* Thicknes ratio function */
  double *BD;                /* Ratio of chord to diameter function */
  double *localChord;        /* Local chord */
  double *DeltaBeta;         /* Twist angle variation function */
  double *sigma;             /* local solidity */
  double *Cl;                /* Lift coefficient function */
  double *Cd;                 /* Drag coefficient function */
  double *manufLambda;       /* Manufacture blade sweep function */
  double *MCA;               /* Midchord alignment function */
}InterpAtZ;

typedef struct{
  double  alt;
  double  T;
  double  p;
  double  c0;
  double  rho;
}AtmData;

typedef struct{
  int NumStations;
  char FN[256];
  d6_t   *BladeVars;
  double *r_R;
  double *Bchord;
  double *Bangle;
  double *dT_dr;
  double *dQ_dr;
  double *AA;                /* Airfoil cross-sectional area */
  double MCA;
}DebugDev;

typedef struct{
  double  D;
  int     NB;
  double  RPM;
  double  Mach_x;
  double  Mach_h;
  double  Mach_t;            /* Tip rotational Mach number */
  double  *Mach_r;           /* sectional relative Mach number */
  double  Mach_tr;           /* Blade tip rotational mach number -MT */
  double  alpha;             /* pitch angle of propeller shaft axis relative to flight direction, rad */
  double  CT;
  double  CP;
  double  Cl;
  double  Cd;
  double  J;
  double  eta;
  double  omega;
  double  Thrust;
  double  Torque;
  double  Re;
  double  Rt;               /* blade tip radius */
  double  bladeAngle;
  double  bladeChord;
  int     HNum;
  BladeGeomData bladeGeom;
  BladeAeroData bladeAero;
  double    T;                /* rotation period */
  double    *t;               /* time */
  double    *Z;               /* Normalized blade radial station - Z [m?] */
  InterpAtZ atZ;
  AtmData   *atm;
  DebugDev  DbDev;
}PropellerData;

typedef struct{
  int numMics;
  int numHarmo;
  double *theta;
  double *y;           /* listener vertical positions */
  double *r;            /* distance observer-propeller */
  double complex **pmL;
  double complex **pmT;
  double **OpT;
  double **OpL;
  double **SPLL;
  double **SPLT;
  double **SPL;

}ListenerMics;


typedef struct{
  int     thetaNum;
  int     phiNum;
  int     NFNum;
  int     FFNum;
  double  *theta;
  double  *theta_r;         /* observer angle in the retarted time reference frame */
  double  *theta_rp;        /* observer angle in the retarted time reference frame */
  double  *phi;             /* tangential angle, arctan(z=y), rad*/
  double  *phi_p;           /* tangential angle, arctan(z=y) relative propeller shaft axis, rad */
  d3_t    *NFcoords;        /* near field */
  d3_t    *FFcoords;        /* far field */
  int     ifRead;
  char    coordsFN[256];
  double  *Y;
  double  *S0;
  double  *S;
  double  *S_r;
  double  **kx;
  double  **phi_s;
  ListenerMics Mics;
}ObserverData;

typedef struct{
  int kx_size;
  double *kx_plt;
  double complex *parabolic_plt;
  double complex *NACA00XX_plt;
  double complex *uniform_plt;
}srcPlot;

typedef struct{
  int plt_size;
  double *kx_plt;
  srcPlot PsiV;
  srcPlot PsiL;
  srcPlot PsiD;
}srcType;




  
  
  
enum runtype{
  GarrickWatkins=0,\
  GutinDeming=1,\
  BarryMagliozzi=2,\
  Hanson=3 
};

//double complex ii=0.0+1.0*I;





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
  /* Debug&Dev alloc-dealloc functions */
  PropellerData allocDbDevInput(PropellerData prop);
  void deallocDbDevInput(PropellerData *prop);
  /* Blade Geom alloc-dealloca functions */
  PropellerData allocBladeGeom(PropellerData prop, d5_t *bladeGeom,int N);
  void deallocBladeGeom(PropellerData prop,int N);

  PropellerData allocBladeAeroData(PropellerData prop, d6_t *bladeAero, int N);
  void deallocBladeAeroData(PropellerData prop, int N);

  void plotDistributions(PropellerData prop,int N);

  /*-------------------------------------------------------------------
              Hanson's method distribution functions                   
                                                                       
  f_D(x) normilized chordwise blade drag function
  f_L(x) normizlied chordwise blade loading function
  H(x)   normizlied chordwise blade thickness distribution
  X      normilized chordwise coordniate
  *------------------------------------------------------------------*/
  /* Set of Normalized thickness Distribution functions H(x) */
  ObserverData getRetard(PropellerData prop, ObserverData obsrvr);
  void deallocObsrvrs(PropellerData prop,ObserverData obsrvr);

  double PsiV(double kx, int m);
  double PsiL(double kx, int m);

  void getNoise(caseData caseIN, PropellerData prop, AtmData atm, ObserverData obsrvr);


  double Hx_parabolic(double X);
  double Hx_NACA00XX(double X);
  double Hx_NACA16seri(double X);
  double fHx(double X);

  /* Lift distribution function - f_L(x) */
  double fLx_uniform(double X);/* unifrom Lift distribution */
  double fLx_parabolic(double X);/* parabolic Lift distribution */
  double fLx(double X);/* Lift distrubution function that chooses a function to simulate */

  /* Drag distribution function - f_D(x) */
  double fDx_uniform(double X);/* unifrom Lift distribution */
  double fDx_parabolic(double X);/* parabolic Lift distribution */
  double fDx(double X);/* Lift distrubution function that chooses a function to simulate */


  double complex cfexp(double x,double kxx);

  double IntPsiV_parabolic(double x);
  double IntPsiV_NACA00XX(double x);
  double IntPsiL_uniform(double x);
  double IntPsiL_parabolic(double x);
  double IntPsiD_uniform(double x);
  double IntPsiD_parabolic(double x);

  /*allocate srce plot */
  srcType allocSourcePlot(srcType src);
  void deallocSourcePlot(srcType src);

  /* Listener/receiver functions */



  /* Mem functions for user defined data types */
  d2_t *d2_tvector(long nl,long nh);
  void free_d2_tvector(d2_t *inputAr,long nl,long nh);

  d3_t *d3_tvector(long nl,long nh);
  void free_d3_tvector(d3_t *inputAr,long nl,long nh);
/*--------------------------- File reader utility functions ----------*/
  d3_t *freader3col(char* fn,int* n);

  d4_t *polarDataReader(char* fn,int* n);

  d5_t *freader5col(char* fn,int* n);

  d6_t *freader6col(char* fn,int* n);

  Airfoil_t *AirfoilDBreader(char* fn,int* n);


/* --------------------------- Math utility functions ----------------*/
  /* Arithmetics */

  /*Interpolation */
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

  void ratint(double xa[],\
              double ya[],\
              int n,\
              double x,\
              double *y);

  void polint(double xa[],\
              double ya[],\
              int n,\
              double x,\
              double *y);

  double getMax(double arr[], int n);

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

/*---------------------------------------------------------------------
                       INTEGRATION FUNCTIONS                           
---------------------------------------------------------------------*/
  double complex qsimp(double (*func)(double), double a, double b, double kxx);
  double complex trapzd(double (*func)(double), double a, double b, double kxx, int n);

  double qgaus(double (*func)(double), double a, double b);
  double complex cqgaus(double (*func)(double), double a, double b,double kxx);

  double integTrapz(double *x, int N);

  double complex integTrapzc(double complex *x, int N);

  double factorial(int n);

  double besselj(int n,double x);

  double Psi_L(double kx);
  double kx(PropellerData prop,\
            double c_R,\
            double r_R,\
            double theta_r,\
            int m);

  double complex **prmsH(PropellerData prop,\
                         AtmData atm,\
                         ObserverData obsrvr,\
                         int m);

  double get_cpu_time(clock_t start,\
                      clock_t end);


