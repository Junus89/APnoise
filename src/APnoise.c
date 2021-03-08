#include "APnoise.h"

/* utility functions */
  void APnoise_error(char error_text[],\
                     const char* routine)
  /*-------------------------------------------------------------------
  Purpose:
          reports error for the APnoise Solver implementation.

  Written by:
             Yunusi Fuerkaiti
             email: y.fuerkaiti@tudelft.nl
  date:
       12.06.2020
  -------------------------------------------------------------------*/
  {
    fprintf(stderr,"                                    \n");
    fprintf(stderr,"APnoise run-time error in routine:\t%s\n",\
      routine);
    fprintf(stderr,"%s\n",error_text);
    fprintf(stderr,"Now exiting to system ...\n");
    exit(1);
  }

  void APnoise_warning(char error_text[],\
                       const char* routine)
  /*-------------------------------------------------------------------
  Purpose:
          reports warning for the APnoise Solver implementation.

  Written by:
             Yunusi Fuerkaiti
             email: y.fuerkaiti@tudelft.nl
  date:
       12.06.2020
  -------------------------------------------------------------------*/
  {
    fprintf(stderr,"                                    \n");
    fprintf(stderr,"APnoise run-time warning in routine:\t%s\n",\
      routine);
    fprintf(stderr,"%s\n",error_text);
    fprintf(stderr,"Please check your routine and run again! \n");
    //exit(1);
  }


/* memory unitilty functions */
  double *dvector(long nl,long nh)
  /*-------------------------------------------------------------------
  Purpose:
          allocates an array of double type.

  Written by:
             Yunusi Fuerkaiti
             email: y.fuerkaiti@tudelft.nl
  date:
       12.06.2020
  -------------------------------------------------------------------*/
  {
    double *outputAr=NULL;
    const char* thisroutine="double *dvector(...)";
    outputAr = (double *)malloc((nh+1)*sizeof(double));
    if(outputAr == NULL)
    {
      APnoise_error("No enough memory space.",thisroutine);
    }
    return outputAr;
  }


  void free_dvector(double *inputAr,\
                    long nl,\
                    long nh)
  {
      free(inputAr);
      return;
  }

  /* complex dvector */
  double complex *dcvector(long nl,long nh)
  /*-------------------------------------------------------------------
  Purpose:
          allocates an array of double complex type.

  Written by:
             Yunusi Fuerkaiti
             email: y.fuerkaiti@tudelft.nl
  date:
       12.06.2020
  -------------------------------------------------------------------*/
  {
    const char* thisroutine="complex *dvector(...)";

    double complex *outputAr=NULL;

    outputAr = (double complex*)malloc((nh+1)*sizeof(double));
    if (outputAr == NULL)
    {
      APnoise_error("No enough memory space.",thisroutine);
    }
    return outputAr;
  }


  void free_dcvector(double complex *inputAr,\
                     long nl,\
                     long nh)
  {
    free(inputAr);
    return;
  }
  /*-----------------------------2D----------------------------------*/

  double **dmatrix(long nrl,\
                   long nrh,\
                   long ncl,\
                   long nch)
  /*-------------------------------------------------------------------
  Purpose:
          allocates a matrix of double type.

  Written by:
             Yunusi Fuerkaiti
             email: y.fuerkaiti@tudelft.nl
  date:
       12.06.2020
  -------------------------------------------------------------------*/
  {
    const char* thisroutine="double **dmatrix(...)";

    int i;
    double **outputAr=NULL;
    outputAr = (double **)malloc((nrh+1)*sizeof(double *));
    if (outputAr == NULL)
    {
      APnoise_error("No enough memory space.",thisroutine);
    }

    for (i = 0; i < nrh+1; i++){
      outputAr[i] = (double *)malloc((nch+1)*sizeof(double));
      if (outputAr[i] == NULL){
      APnoise_error("No enough memory space.",thisroutine);
      }
    }
    return outputAr;
  }
  /* */
  void free_dmatrix(double **inputAr,\
                    long nrl,\
                    long nrh,\
                    long ncl,\
                    long nch)
  /*      
   */
  {
      int i;
      for (i = 0; i <= nrh; i++) free(inputAr[i]);
      free(inputAr);              
      return;
  }

/* */
/* complex dmatrix */
/* complex dmatrix */
  double complex **dcmatrix(long nrl,\
                            long nrh,\
                            long ncl,\
                            long nch)
  /*-------------------------------------------------------------------
  Purpose:
          allocates a matrix of double type.

  Written by:
             Yunusi Fuerkaiti
             email: y.fuerkaiti@tudelft.nl
  date:
       12.06.2020
  -------------------------------------------------------------------*/
  {
    const char* thisroutine="double complex **dcmatrix(...)";
    int i;
    double complex **outputAr=NULL;
    outputAr = (double complex**)malloc((nrh+1)*sizeof(double complex*));
    if (outputAr == NULL)
    {
    APnoise_error("No enough memory space.",thisroutine);
    }

    for (i = 0; i < nrh+1; i++){
      outputAr[i] = (double complex*)malloc((nch+1)*sizeof(double complex));
      if (outputAr[i] == NULL){
        APnoise_error("No enough memory space.",thisroutine);
      }
    }
    return outputAr;
  }
/* */

  /* */
  void free_dcmatrix(double complex **inputAr,\
                     long nrl,\
                     long nrh,\
                     long ncl,\
                     long nch)
  /*      
   */
  {
      int i;
      for (i = 0; i <= nrh; i++)  free(inputAr[i]);
      free(inputAr);              
      return;
  }

  /* */
  d2_t *d2_tvector(long nl,long nh)
  /*-------------------------------------------------------------------
  Purpose:
          allocate a vector of d2_t type.

  Written by:
             Yunusi Fuerkaiti
             email: y.fuerkaiti@tudelft.nl
  date:
       12.06.2020
  -------------------------------------------------------------------*/
  {
    const char* thisroutine="double complex ***dc3dtensor(...)";
    d2_t *outputAr=NULL;

    outputAr=(d2_t *)malloc((nh+1)*sizeof(d2_t));
    if (outputAr==NULL){
      APnoise_error("No enough memory space.",thisroutine);
    }
    return outputAr;
  }


  void free_d2_tvector(d2_t *inputAr,\
                       long nl,\
                       long nh)
  {
      free(inputAr);
      return;
  }

  d2_t **d2_tmatrix(long nrl,\
                    long nrh,\
                    long ncl,\
                    long nch)
    /*-------------------------------------------------------------------
    Purpose:
            allocate a matrix of d2_t type.

    Written by:
               Yunusi Fuerkaiti
               email: y.fuerkaiti@tudelft.nl
    date:
         12.06.2020
    -------------------------------------------------------------------*/
  {
    int i;
    const char* thisroutine="d2_t **d2_tmatrix(...)";
    d2_t **outputAr=NULL;
    outputAr=(d2_t **)malloc((nrh+1)*sizeof(d2_t *));
    if (outputAr==NULL)
    {
      APnoise_error("No enough memory space.",thisroutine);
    }

    for (i=0;i<nrh+1;i++){
      outputAr[i]=(d2_t *)malloc((nch+1)*sizeof(d2_t));
      if (outputAr[i] == NULL){
        APnoise_error("No enough memory space.",thisroutine);
      }
    }
    return outputAr;
  }

  /* */
  void free_d2_tmatrix(d2_t **inputAr,\
                       long nrl,\
                       long nrh,\
                       long ncl,\
                       long nch)
  /*      
   */
  {
    int i;
    for (i=0;i<=nrh;i++)  free(inputAr[i]);
    free(inputAr);
    return;
  }


  double ***d3dtensor(long nrl,
                      long nrh,\
                      long ncl,\
                      long nch,\
                      long ndl,\
                      long ndh)
  /*-------------------------------------------------------------------
  Purpose:
          allocates a tensor(3D vector) of double type.

  Written by:
             Yunusi Fuerkaiti
             email: y.fuerkaiti@tudelft.nl
  date:
       12.06.2020
  -------------------------------------------------------------------*/
  {
    const char* thisroutine="double ***d3dtensor(...)";
    int i,j;
    double ***outputAr=NULL;
    outputAr = (double ***)malloc((nrh+1)*sizeof(double **));
    if (outputAr == NULL)
    {
        APnoise_error("No enough memory space.",thisroutine);
        exit(1);
    }

    for (i = 0; i < nrh+1; i++){
      outputAr[i] = (double **)malloc((nch+1)*sizeof(double*));
      if (outputAr[i] == NULL)
              {
      APnoise_error("No enough memory space.",thisroutine);
              }
      }

    for (i = 0; i<nrh+1;i++){
      for (j = 0;j<nch+1;j++){
        outputAr[i][j]=(double *)malloc((ndh+1)*sizeof(double));
        if (outputAr[i][j]==NULL){
          APnoise_error("No enough memory space.",thisroutine);
        }
      }
    }
    return outputAr;
  }

  void free_d3dtensor(double ***inputAr,
                      long nrl,\
                      long nrh,\
                      long ncl,\
                      long nch,\
                      long ndl,\
                      long ndh)
  /*      
   */
  {
    int i,j;

    for (i = 0; i < nrh+1; i++) {
      for (j = 0; j < nch+1; j++){
        free(inputAr[i][j]);
      }
    }
    for (i = 0; i < nrh+1; i++) free(inputAr[i]);
    free(inputAr);
    return;
  }

  d3_t *d3_tvector(long nl,long nh)
    /*-------------------------------------------------------------------
    Purpose:
            allocate a vector of d3_t type.

    Written by:
               Yunusi Fuerkaiti
               email: y.fuerkaiti@tudelft.nl
    date:
         12.06.2020
    -------------------------------------------------------------------*/
  {
    const char* thisroutine="d3_t **d3_tvector(...)";
    d3_t *outputAr=NULL; 
    outputAr = (d3_t *)malloc((nh+1)*sizeof(d3_t));
    if (outputAr == NULL)
    {
      APnoise_error("No enough memory space.",thisroutine);
    }
    return outputAr;
  }


  void free_d3_tvector(d3_t *inputAr,long nl,long nh)
  {
      free(inputAr);
      return;
  }

  PropellerData allocDbDevInput(PropellerData prop)
  {
  /*
  */
  prop.DbDev.MCA=0.00; /* randomly set */

  prop.DbDev.BladeVars=NULL; prop.DbDev.BladeVars=freader6col(prop.DbDev.FN,&prop.DbDev.NumStations);/* free later */
  prop.DbDev.r_R=NULL; prop.DbDev.r_R=dvector(0,prop.DbDev.NumStations);
  prop.DbDev.Bchord=NULL; prop.DbDev.Bchord=dvector(0,prop.DbDev.NumStations);
  prop.DbDev.Bangle=NULL; prop.DbDev.Bangle=dvector(0,prop.DbDev.NumStations);
  prop.DbDev.AA=NULL; prop.DbDev.AA=dvector(0,prop.DbDev.NumStations);
  prop.DbDev.dT_dr=NULL; prop.DbDev.dT_dr=dvector(0,prop.DbDev.NumStations);
  prop.DbDev.dQ_dr=NULL; prop.DbDev.dQ_dr=dvector(0,prop.DbDev.NumStations);
/* */
  prop.Mach_r=NULL;prop.Mach_r=dvector(0,prop.DbDev.NumStations);


  for(int i=0;i<prop.DbDev.NumStations;i++){
   prop.DbDev.r_R[i]=prop.DbDev.BladeVars[i].x[0];
   prop.DbDev.Bchord[i]=prop.DbDev.BladeVars[i].x[1];
   prop.DbDev.Bangle[i]=prop.DbDev.BladeVars[i].x[2];
   prop.DbDev.AA[i]=prop.DbDev.BladeVars[i].x[3]*0.6853;
   prop.DbDev.dT_dr[i]=prop.DbDev.BladeVars[i].x[4];
   prop.DbDev.dQ_dr[i]=prop.DbDev.BladeVars[i].x[5];
   /* 
   */
   prop.Mach_r[i]=sqrt(pow(prop.Mach_x,2)+pow(prop.DbDev.r_R[i],2)*pow(prop.Mach_t,2));

 }

  if(DEBUG==0){
     printf("Dev FN : %s\n",prop.DbDev.FN);
     for(int i=0;i<prop.DbDev.NumStations;i++){
       printf("AA[%d] = %lf\n",i,prop.DbDev.AA[i]);
     }
  }
  return prop;
  }

  void deallocDbDevInput(PropellerData *prop)
  {
  const char* thisroutine="deallocDbDevInput";
  int code=0;
  /*
  */
  free(prop->DbDev.BladeVars);//prop->DbDev.BladeVars=NULL;
  free_dvector(prop->DbDev.r_R,0,prop->DbDev.NumStations);//prop->DbDev.r_R=NULL;
  free_dvector(prop->DbDev.Bchord,0,prop->DbDev.NumStations);
  free_dvector(prop->DbDev.Bangle,0,prop->DbDev.NumStations);
  free_dvector(prop->DbDev.AA,0,prop->DbDev.NumStations);
  free_dvector(prop->DbDev.dT_dr,0,prop->DbDev.NumStations);
  free_dvector(prop->DbDev.dQ_dr,0,prop->DbDev.NumStations);
  free_dvector(prop->Mach_r,0,prop->DbDev.NumStations);


  if(code!=0)APnoise_error(" error in routine %s\n",thisroutine);
  return;
  }



  PropellerData allocBladeGeom(PropellerData prop, d5_t *bladeGeom,int N)
  {
  /*
  */
  prop.bladeGeom.r_R=NULL;prop.bladeGeom.r_R=dvector(0,prop.bladeGeom.NumStations);
  prop.bladeGeom.r=NULL;prop.bladeGeom.r=dvector(0,prop.bladeGeom.NumStations);
  prop.bladeGeom.c_R=NULL;prop.bladeGeom.c_R=dvector(0,prop.bladeGeom.NumStations);
  prop.bladeGeom.c=NULL;prop.bladeGeom.c=dvector(0,prop.bladeGeom.NumStations);
  prop.bladeGeom.twist=NULL;prop.bladeGeom.twist=dvector(0,prop.bladeGeom.NumStations);
  prop.bladeGeom.thickns=NULL;prop.bladeGeom.thickns=dvector(0,prop.bladeGeom.NumStations);
  prop.bladeGeom.c_D=NULL;prop.bladeGeom.c_D=dvector(0,prop.bladeGeom.NumStations);
  prop.bladeGeom.lambdaManufac=NULL;prop.bladeGeom.lambdaManufac=dvector(0,prop.bladeGeom.NumStations);


  for(int i=0;i<prop.bladeGeom.NumStations;i++){
    prop.bladeGeom.r_R[i]=bladeGeom[i].x[0];
    prop.bladeGeom.r[i]=bladeGeom[i].x[1];
    prop.bladeGeom.c_R[i]=bladeGeom[i].x[2];
    prop.bladeGeom.c[i]=bladeGeom[i].x[3];
    prop.bladeGeom.twist[i]=bladeGeom[i].x[4];
    prop.bladeGeom.thickns[i]=prop.bladeGeom.c_R[i]*(prop.D/2)*prop.bladeGeom.thicknsRatio;
    prop.bladeGeom.c_D[i]=prop.bladeGeom.c_R[i]*(prop.D/2)/prop.D;
    prop.bladeGeom.lambdaManufac[i]=prop.bladeGeom.r_R[i]*0.0;
    if(DEBUG==0){
      //printf("bladeGeom.r_R[%d] = %lf, prop.bladeGeom.c_R[%d] = %lf\n",i,prop.bladeGeom.r_R[i],i,prop.bladeGeom.c_R[i]);
      printf("bladeGeom.r_R[%d] = %lf, prop.bladeGeom.thickns[%d] = %lf, thicknsRati=%lf\n",i,prop.bladeGeom.r_R[i],i,prop.bladeGeom.thickns[i]\
      ,prop.bladeGeom.thicknsRatio);
    }
  }
  prop.Mach_tr=sqrt(prop.Mach_t*prop.Mach_t-prop.Mach_x*prop.Mach_x);
  printf("prop.Mach_tr = %lf\n",prop.Mach_tr);
  prop.T=2*PI/prop.omega/prop.NB;
  int d=500;
  prop.t=NULL; prop.t=linspace(-prop.T,prop.T,d);
  prop.Rt=prop.D/2;
  //int N=100;
  prop.bladeGeom.hub_r=0.0625; /*needs to be read from the input file */
  double RH=prop.bladeGeom.hub_r*2;/*hub radius */
  double *rR=NULL;rR=linspace(RH,prop.Rt,N);
  prop.Z=NULL;prop.Z=dvector(0,N);
  prop.Mach_r=NULL;prop.Mach_r=dvector(0,N);
  /* allocate intereplation functions */

  prop.atZ.tb=NULL; prop.atZ.tb=dvector(0,N); d3_t tbINT;
  prop.atZ.BD=NULL; prop.atZ.BD=dvector(0,N); d3_t BDINT;
  prop.atZ.localChord=NULL; prop.atZ.localChord=dvector(0,N);
  prop.atZ.DeltaBeta=NULL; prop.atZ.DeltaBeta=dvector(0,N); d3_t DeltaBetaINT;
  prop.atZ.sigma=NULL; prop.atZ.sigma=dvector(0,N);
  prop.atZ.manufLambda=NULL; prop.atZ.manufLambda=dvector(0,N); d3_t manufLambdaINT;
  prop.atZ.MCA=NULL; prop.atZ.MCA=dvector(0,N); d3_t MCAINT;

  for(int i=0;i<N;i++){
    prop.Z[i]=rR[i]/prop.Rt;
    prop.Mach_r[i]=sqrt(pow(prop.Mach_x,2)+pow(prop.Z[i],2)*pow(prop.Mach_t,2));

    ratint(prop.bladeGeom.r_R,\
           prop.bladeGeom.thickns,\
           prop.bladeGeom.NumStations,\
           prop.Z[i],\
           &prop.atZ.tb[i]); /* tb */

    ratint(prop.bladeGeom.r_R,\
           prop.bladeGeom.c_D,\
           prop.bladeGeom.NumStations,\
           prop.Z[i],\
           &prop.atZ.BD[i]); /* BD */
    prop.atZ.localChord[i]=prop.atZ.BD[i]*prop.D;
    prop.atZ.sigma[i]=prop.NB*prop.atZ.localChord[i]/(2*PI*rR[i]);
    ratint(prop.bladeGeom.r_R,\
           prop.bladeGeom.twist,\
           prop.bladeGeom.NumStations,\
           prop.Z[i],\
           &prop.atZ.DeltaBeta[i]); /* DeltaBeta*/
    ratint(prop.bladeGeom.r_R,\
           prop.bladeGeom.lambdaManufac,\
           prop.bladeGeom.NumStations,\
           prop.Z[i],\
           &prop.atZ.manufLambda[i]); /* manufLambda */
    prop.atZ.MCA[i]=tan((prop.atZ.manufLambda[i]*PI/180.0))*prop.Z[i];
  }
  if(DEBUG==0){
    for(int i=0;i<N;i++){
     // printf("Numstations  = %d\n",prop.bladeGeom.NumStations);
    //printf("prop.Z[%d] = %lf\n",i,prop.Z[i]);
    //printf("R[%d]=%lf\n",i,rR[i]);
    //printf("Mach_r[%d]=%lf\n",i,prop.Mach_r[i]);
    //printf("Mach_x=%lf, Mach_t = %lf\n",prop.Mach_x,prop.Mach_t);
    //printf("tbINT[%d] = %lf\n",i,prop.atZ.BD[i]);
     printf("Z[%d] = %lf thikns = %lf tb = %lf\n",i,prop.Z[i],prop.bladeGeom.thickns[i],prop.atZ.tb[i]);
  }
  }
  return prop;
  }

  void deallocBladeGeom(PropellerData prop,int N)
  {
  free_dvector(prop.bladeGeom.r_R,0,prop.bladeGeom.NumStations);
  free_dvector(prop.bladeGeom.r,0,prop.bladeGeom.NumStations);
  free_dvector(prop.bladeGeom.c_R,0,prop.bladeGeom.NumStations);
  free_dvector(prop.bladeGeom.c,0,prop.bladeGeom.NumStations);
  free_dvector(prop.bladeGeom.twist,0,prop.bladeGeom.NumStations);
  free_dvector(prop.bladeGeom.thickns,0,prop.bladeGeom.NumStations);
  free_dvector(prop.bladeGeom.c_D,0,prop.bladeGeom.NumStations);
  free_dvector(prop.bladeGeom.lambdaManufac,0,prop.bladeGeom.NumStations);
  free_dvector(prop.Z,0,N);
  free_dvector(prop.Mach_r,0,N);

  free_dvector(prop.atZ.tb,0,N);
  free_dvector(prop.atZ.BD,0,N);
  free_dvector(prop.atZ.DeltaBeta,0,N);
  free_dvector(prop.atZ.manufLambda,0,N);
  free_dvector(prop.atZ.MCA,0,N);


  }

  PropellerData allocBladeAeroData(PropellerData prop, d6_t *bladeAero, int N)
  {
  /*
  */
  prop.bladeAero.Cl=NULL;prop.bladeAero.Cl=dvector(0,prop.bladeAero.NumStations);
  prop.bladeAero.Cd=NULL;prop.bladeAero.Cd=dvector(0,prop.bladeAero.NumStations);
  prop.bladeAero.AoA=NULL;prop.bladeAero.AoA=dvector(0,prop.bladeAero.NumStations);
  prop.bladeAero.Veff=NULL;prop.bladeAero.Veff=dvector(0,prop.bladeAero.NumStations);
  prop.bladeAero.Ct=NULL;prop.bladeAero.Ct=dvector(0,prop.bladeAero.NumStations);
  prop.bladeAero.Cq=NULL;prop.bladeAero.Cq=dvector(0,prop.bladeAero.NumStations);

  for(int i=0;i<prop.bladeAero.NumStations;i++){
    prop.bladeAero.Cl[i]=bladeAero[i].x[0];
    prop.bladeAero.Cd[i]=bladeAero[i].x[1];
    prop.bladeAero.AoA[i]=bladeAero[i].x[2];
    prop.bladeAero.Veff[i]=bladeAero[i].x[3];
    prop.bladeAero.Ct[i]=bladeAero[i].x[4];
    prop.bladeAero.Cq[i]=bladeAero[i].x[5];

  }

  prop.atZ.Cl=NULL; prop.atZ.Cl=dvector(0,N); d3_t ClINT;
  prop.atZ.Cd=NULL; prop.atZ.Cd=dvector(0,N); d3_t CdINT;

  if(DEBUG==0){
    printf("numstations = %d\n",prop.bladeAero.NumStations);
    for(int k=0;k<prop.bladeAero.NumStations;k++) printf("r_R[%d] = %lf, prop.bladeAero.cl[%d] = %lf\n",\
            k,prop.bladeGeom.r_R[k],k,prop.bladeAero.Cl[k]);
  }
  for(int i=0;i<N;i++){
    ratint(prop.bladeGeom.r_R,\
           prop.bladeAero.Cl,\
           prop.bladeAero.NumStations,\
           prop.Z[i],\
           &prop.atZ.Cl[i]);
    ratint(prop.bladeGeom.r_R,\
           prop.bladeAero.Cd,\
           prop.bladeAero.NumStations,\
           prop.Z[i],\
           &prop.atZ.Cd[i]);
  
  }
  if(DEBUG==0){
     //for(int i=0;i<N;i++) printf("Z[%d] = %lf Cl = %lf Cd = %lf\n",i,prop.Z[i],prop.atZ.Cl[i],prop.atZ.Cd[i]);
     //for(int i=0;i<40;i++) printf("Z[%d] = %lf thicknss = %lf tb = %lf\n",i,prop.Z[i],prop.atZ.tb[i]);

  }


  return prop;
  }

  void deallocBladeAeroData(PropellerData prop, int N)
  {
  free_dvector(prop.bladeAero.Cl,0,prop.bladeAero.NumStations);
  free_dvector(prop.bladeAero.Cd,0,prop.bladeAero.NumStations);
  free_dvector(prop.bladeAero.AoA,0,prop.bladeAero.NumStations);
  free_dvector(prop.bladeAero.Veff,0,prop.bladeAero.NumStations);
  free_dvector(prop.bladeAero.Ct,0,prop.bladeAero.NumStations);
  free_dvector(prop.bladeAero.Cq,0,prop.bladeAero.NumStations);

  free_dvector(prop.atZ.Cl,0,N);
  free_dvector(prop.atZ.Cd,0,N);

  }

  void plotDistributions(PropellerData prop,int N)
  {
  /* */
  FILE *Ftb,*FBD,*FDeltaBeta,*FCl,*FCd;
  Ftb=fopen("Thickns_distrib.txt","w");
  FBD=fopen("BladeWidth_distrib.txt","w");
  FDeltaBeta=fopen("Twist_distrib.txt","w");
  FCl=fopen("Cl_distrib.txt","w");
  FCd=fopen("Cd_distrib.txt","w");

  fprintf(Ftb,"#r_R\t thickns\n");
  fprintf(FBD,"#r_R\t Blade width\n");
  fprintf(FDeltaBeta,"#r_R\t twist\n");
  fprintf(FCl,"#r_R\t Cl\n");
  fprintf(FCd,"#r_R\t Cd\n");

  for(int i=0;i<N;i++){
    fprintf(Ftb,"%4.8f %4.8f\n",prop.Z[i],prop.atZ.tb[i]);
    fprintf(FBD,"%4.8f %4.8f\n",prop.Z[i],prop.atZ.BD[i]);
    fprintf(FDeltaBeta,"%4.8f %4.8f\n",prop.Z[i],prop.atZ.DeltaBeta[i]);
    fprintf(FCl,"%4.8f %4.8f\n",prop.Z[i],prop.atZ.Cl[i]);
    fprintf(FCd,"%4.8f %4.8f\n",prop.Z[i],prop.atZ.Cd[i]);
  }
  fclose(Ftb),fclose(FBD),fclose(FDeltaBeta),fclose(FCl),fclose(FCd);
  
  }

  /*-------------------------------------------------------------------
                                                                       
              Hanson's method distribution functions                   
                                                                       
  f_D(x) normilized chordwise blade drag function
  f_L(x) normizlied chordwise blade loading function
  H(x)   normizlied chordwise blade thickness distribution
  X      normilized chordwise coordniate
                                                                       
  *------------------------------------------------------------------*/

  ObserverData getRetard(PropellerData prop, ObserverData obsrvr)
  {
  /*
  */
  const char* thisroutine="getRetard";

  obsrvr.theta_r=NULL;obsrvr.theta_r=dvector(0,obsrvr.FFNum);
  obsrvr.theta_rp=NULL;obsrvr.theta_rp=dvector(0,obsrvr.FFNum);
  obsrvr.S_r=NULL;obsrvr.S_r=dvector(0,obsrvr.FFNum);
  obsrvr.phi=NULL;obsrvr.phi=dvector(0,obsrvr.FFNum);
  obsrvr.phi_p=NULL;obsrvr.phi_p=dvector(0,obsrvr.FFNum);

  obsrvr.kx=NULL;obsrvr.kx=dmatrix(0,prop.DbDev.NumStations,0,obsrvr.FFNum);
  obsrvr.phi_s=NULL;obsrvr.phi_s=dmatrix(0,prop.DbDev.NumStations,0,obsrvr.FFNum);

  double theta_rpm=0.0, temp=0.0, tempPprm=0.0;

  for(int i=0;i<obsrvr.FFNum;i++){
    obsrvr.phi[i]=atan2(obsrvr.FFcoords[i].x[2],obsrvr.FFcoords[i].x[1]);
    temp=cos(obsrvr.theta[i])*sqrt(1-pow(prop.Mach_x,2)*pow(sin(obsrvr.theta[i]),2)) \
        +prop.Mach_x*pow(sin(obsrvr.theta[i]),2);
    obsrvr.theta_r[i]=acos(temp);
    obsrvr.S_r[i]=obsrvr.Y[i]/sin(obsrvr.theta_r[i]);
    theta_rpm=cos(obsrvr.theta_r[i])*cos(prop.alpha)+sin(obsrvr.theta_r[i]) \
             *sin(obsrvr.phi[i])*sin(prop.alpha);
    obsrvr.theta_rp[i]=acos(theta_rpm);
    tempPprm=sin(obsrvr.theta_r[i])/sin(obsrvr.theta_rp[i])*cos(obsrvr.phi[i]);
    obsrvr.phi_p[i]=acos(tempPprm);
    /* kx and phi_s */
    for(int j=0;j<prop.DbDev.NumStations;j++){
      /* Note this kx must be multiplied with m in the main routine */
      obsrvr.kx[j][i]=2*prop.NB*prop.DbDev.Bchord[j]*prop.Mach_t \
        /(prop.Mach_r[j]*(1-prop.Mach_x*cos(obsrvr.theta_r[i])));
      obsrvr.phi_s[j][i]=2*prop.NB*prop.Mach_t/(prop.Mach_r[j] \
        *(1-prop.Mach_x*cos(obsrvr.theta_r[i])))*prop.DbDev.MCA/prop.D;
    }
  }

  if(DEBUG==1){
    printf("\n\n Debug info from routine %s\n\n",thisroutine);
    for(int i=0;i<prop.DbDev.NumStations;i++){
      printf("Mach_r[%d] = %lf\n",i,prop.Mach_r[i]);
    }
    printf("obsrvr.theta[0]= %lf, theta_r[0] = %lf\n",obsrvr.theta[0],obsrvr.theta_r[0]);
    printf("Y[0] = %lf, S_r[0] = %lf,S[0]= %lf, S0[0] = %lf\n",obsrvr.Y[0],obsrvr.S_r[0],obsrvr.S[0],obsrvr.S0[0]);
    printf("prop.alpha = %lf\n",prop.alpha);
    printf("\n\n");
    printf("Mach_x = %lf, Mach_t = %lf\n",prop.Mach_x,prop.Mach_t);
    printf("NB = %d, RPM = %lf, MCA=%lf\n",prop.NB,prop.RPM,prop.DbDev.MCA);
  }
  return obsrvr;
  }

  void getNoise(caseData caseIN \
               ,PropellerData prop \
               ,AtmData atm \
               ,ObserverData obsrvr)
  {
  /*
  */

  rmdir("SPL",0777);
  mkdir("SPL",0777);

/*
  FILE *fsplm,*fsplr;
  char fnm[256];strcpy(fnm,caseIN.caseName);
  char fnr[256];strcpy(fnr,caseIN.caseName);

  strcat(fnm,"_SPLHansonm.txt");
  strcat(fnr,"_SPLHansonr.txt");
  fsplm=fopen(fnm,"w");
  fsplr=fopen(fnr,"w");
*/
  double p_ref=2e-5;
  double complex ii=(0.0+1.0*I);
  double freq=prop.NB*prop.omega/2/PI;

  obsrvr.Mics.numHarmo=10;/*this needs to be read from the input file */
  obsrvr.Mics.numMics=obsrvr.FFNum;
  obsrvr.Mics.pmL=NULL; obsrvr.Mics.pmL=dcmatrix(0,obsrvr.Mics.numHarmo,0,obsrvr.FFNum);
  obsrvr.Mics.pmT=NULL; obsrvr.Mics.pmT=dcmatrix(0,obsrvr.Mics.numHarmo,0,obsrvr.FFNum);
  obsrvr.Mics.OpL=NULL; obsrvr.Mics.OpL=dmatrix(0,obsrvr.Mics.numHarmo,0,obsrvr.FFNum);
  obsrvr.Mics.OpT=NULL; obsrvr.Mics.OpT=dmatrix(0,obsrvr.Mics.numHarmo,0,obsrvr.FFNum);
  obsrvr.Mics.SPLL=NULL; obsrvr.Mics.SPLL=dmatrix(0,obsrvr.Mics.numHarmo,0,obsrvr.FFNum);
  obsrvr.Mics.SPLT=NULL; obsrvr.Mics.SPLT=dmatrix(0,obsrvr.Mics.numHarmo,0,obsrvr.FFNum);
  obsrvr.Mics.SPL=NULL; obsrvr.Mics.SPL=dmatrix(0,obsrvr.Mics.numHarmo,0,obsrvr.FFNum);

  double complex integRoot2TipL=(0.0+0.0*I),integRoot2TipT=(0.0+0.0*I);
  double complex *integRTL=NULL; integRTL=dcvector(0,prop.DbDev.NumStations);
  double complex *integRTT=NULL; integRTT=dcvector(0,prop.DbDev.NumStations);

  double PsiLm=0.0, PsiVm=0.0;
  for(int k=0;k<prop.DbDev.NumStations;k++){
    prop.DbDev.Bangle[k]=prop.DbDev.Bangle[k]*PI/180.0;
  }

  /* Loading term */
  for(int m=1;m<prop.HNum+1;m++){
    for(int j=0;j<obsrvr.FFNum;j++){
      for(int k=0;k<prop.DbDev.NumStations;k++){
        PsiLm=PsiL(obsrvr.kx[j][k],m);
        integRTL[k]=((cos(obsrvr.theta_rp[j])/(1-prop.Mach_x*cos(obsrvr.theta_r[j]))*prop.DbDev.dT_dr[k] \
          -1/(pow(prop.DbDev.r_R[k],2)*prop.Mach_t*prop.D/2)*prop.DbDev.dQ_dr[k])*cexp(ii*obsrvr.phi_s[k][j]) \
          *besselj(m*prop.NB,m*prop.NB*prop.DbDev.r_R[k]*prop.Mach_t*sin(obsrvr.theta_rp[j]) \
          /(1-prop.Mach_x*cos(obsrvr.theta_r[j])))*PsiLm);
      }
      integRoot2TipL=integTrapzc(integRTL,prop.DbDev.NumStations);
      obsrvr.Mics.pmL[m][j]=ii*prop.NB*prop.Mach_t*sin(obsrvr.theta_r[j]) \
        *cexp(ii*m*prop.NB*(prop.omega*obsrvr.S_r[j]/atm.c0+(obsrvr.phi_p[j]-PI/2))) \
        /(2*sqrt(2)*PI*obsrvr.Y[j]*prop.D/2*(1-prop.Mach_x*cos(obsrvr.theta_r[j])))*integRoot2TipL;
      obsrvr.Mics.OpL[m][j]=1*sqrt(pow(creal(obsrvr.Mics.pmL[m][j]),2)+pow(cimag(obsrvr.Mics.pmL[m][j]),2));
      obsrvr.Mics.SPLL[m][j]=20*log10(fabs(obsrvr.Mics.OpL[m][j])/p_ref);
    }
  }
  /* Thickness term */
  for(int m=1;m<prop.HNum+1;m++){
    for(int j=0;j<obsrvr.FFNum;j++){
      for(int k=0;k<prop.DbDev.NumStations;k++){
        PsiVm=PsiV(obsrvr.kx[j][k],m);/* Note AA instead of (h/b) */
        integRTT[k]=(pow(prop.Mach_r[k],2)*(prop.DbDev.AA[k]/1)*cexp(ii*obsrvr.phi_s[k][j]) \
          *besselj(m*prop.NB,m*prop.NB*prop.DbDev.r_R[k]*prop.Mach_t*sin(obsrvr.theta_rp[j]) \
          /(1-prop.Mach_x*cos(obsrvr.theta_r[j])))*pow((obsrvr.kx[j][k]*m),2)*PsiVm);
      }
      integRoot2TipT=integTrapzc(integRTT,prop.DbDev.NumStations);
      obsrvr.Mics.pmT[m][j]=-atm.rho*pow(atm.c0,2)*prop.NB*sin(obsrvr.theta_r[j]) \
        *cexp(ii*m*prop.NB*(prop.omega*obsrvr.S_r[j]/atm.c0+(obsrvr.phi_p[j]-PI/2))) \
        /(4*sqrt(2)*PI*(obsrvr.Y[j]/prop.D)*(1-prop.Mach_x*cos(obsrvr.theta_r[j])))*integRoot2TipT;
      obsrvr.Mics.OpT[m][j]=1*sqrt(pow(creal(obsrvr.Mics.pmT[m][j]),2)+pow(cimag(obsrvr.Mics.pmT[m][j]),2));
      obsrvr.Mics.SPLT[m][j]=20*log10(fabs(obsrvr.Mics.OpT[m][j])/p_ref);
    }
  }
  /* Total */
    for(int j=0;j<obsrvr.FFNum;j++){
      FILE *fp[j], *fpr[j];
      char fn[2048], fnr[2048];
      strcpy(fn,"SPL/"); strcpy(fnr,"SPL/");
      strcat(fn,caseIN.caseName); strcat(fnr,caseIN.caseName);
      sprintf(fn+strlen(fn),"_SPLH_Mic%d.txt",j); sprintf(fnr+strlen(fnr),"_SPLHr_Mic%d.txt",j);
      /* Total */
      fp[j]=fopen(fn,"w");
      fprintf(fp[j],"%5s %5s %5s %s \n","# f [hz] ", "SPL_Th [dB]", "SPL_L [dB]", "SPL [dB]");
      for(int m=1;m<prop.HNum+1;m++){
        obsrvr.Mics.SPL[m][j]=20*log10(fabs(obsrvr.Mics.OpT[m][j]+obsrvr.Mics.OpL[m][j])/p_ref);

/*      if(m==1) fprintf(fsplr,"%4.4f %4.4f %4.4f %4.4f\n",(180-obsrvr.theta[j]*180.0/PI),obsrvr.Mics.SPLT[1][j] \
        ,obsrvr.Mics.SPLL[1][j],obsrvr.Mics.SPL[1][j]);
       fprintf(fsplm,"%4.6f %4.6f %4.6f %4.6f\n",m*freq,obsrvr.Mics.SPLT[m][0],obsrvr.Mics.SPLL[m][0] \
      ,obsrvr.Mics.SPL[m][0]);
*/
       fprintf(fp[j],"%4.6f %4.6f %4.6f %4.6f\n",m*freq,obsrvr.Mics.SPLT[m][j],obsrvr.Mics.SPLL[m][j] \
        ,obsrvr.Mics.SPL[m][j]);
    }
    fclose(fp[j]);
    fpr[j]=fopen(fnr,"w");
    fprintf(fpr[j],"%5s %5s %5s %s \n","Theta [deg]", "SPL_Th [dB]", "SPL_L [dB]", "SPL [dB]");
    for(int m=1;m<prop.HNum+1;m++){
       fprintf(fpr[j],"%4.6f %4.6f %4.6f %4.6f\n",(180-obsrvr.theta[j]*180.0/PI),obsrvr.Mics.SPLT[m][j] \
        ,obsrvr.Mics.SPLL[m][j],obsrvr.Mics.SPL[m][j]);
    }

    fclose(fpr[j]);
  }
/*  fclose(fsplm);
  fclose(fsplr);
*/
  if(DEBUG==1){
    printf("obsrvr numbers = %d, num mics = %d\n",obsrvr.FFNum,obsrvr.Mics.numMics);
    printf("blade angle = %lf\n",prop.DbDev.Bangle[10]);
  }
  free_dcvector(integRTL,0,prop.DbDev.NumStations);
  free_dcvector(integRTT,0,prop.DbDev.NumStations);

  free_dcmatrix(obsrvr.Mics.pmL,0,obsrvr.Mics.numHarmo,0,obsrvr.FFNum);
  free_dcmatrix(obsrvr.Mics.pmT,0,obsrvr.Mics.numHarmo,0,obsrvr.FFNum);
  free_dmatrix(obsrvr.Mics.OpL,0,obsrvr.Mics.numHarmo,0,obsrvr.FFNum);
  free_dmatrix(obsrvr.Mics.OpT,0,obsrvr.Mics.numHarmo,0,obsrvr.FFNum);
  free_dmatrix(obsrvr.Mics.SPLL,0,obsrvr.Mics.numHarmo,0,obsrvr.FFNum);
  free_dmatrix(obsrvr.Mics.SPLT,0,obsrvr.Mics.numHarmo,0,obsrvr.FFNum);
  free_dmatrix(obsrvr.Mics.SPL,0,obsrvr.Mics.numHarmo,0,obsrvr.FFNum);


  return;
  }

  void deallocObsrvrs(PropellerData prop,ObserverData obsrvr)
  {
  /*
  */
  free_dvector(obsrvr.theta_r,0,obsrvr.FFNum);
  free_dvector(obsrvr.theta_rp,0,obsrvr.FFNum);
  free_dvector(obsrvr.S_r,0,obsrvr.FFNum);
  free_dvector(obsrvr.phi,0,obsrvr.FFNum);
  free_dvector(obsrvr.phi_p,0,obsrvr.FFNum);

  free_dmatrix(obsrvr.kx,0,prop.DbDev.NumStations,0,obsrvr.FFNum);
  free_dmatrix(obsrvr.phi_s,0,prop.DbDev.NumStations,0,obsrvr.FFNum);
  return;
  }

  double PsiV(double kx, int m)
  {
  /* This function evaluates the Psiv term in the Hanson's formulation 
  For details please refer to : Kotwicz Herniczek et al., 2017
  
  */
  double psiv;
  double kxm=kx*m; /* this should be done in this way, as the input
  kx is actually is kx_literature/m; */
  if(kxm==0.0){
    psiv=2/3;
  }
  else{
    psiv=8/pow(kxm,2)*(2/kxm*sin(kxm/2)-cos(kxm/2));
  }
  return psiv;
  }

  double PsiL(double kx, int m)
  {
  /* This function evaluates the PsiL term in the Hanson's formulation 
  For details please refer to : Kotwicz Herniczek et al., 2017
  */
  double psil;
  double kxm=kx*m;
  if(kxm==0.0){
    psil=1.0;
  }
  else{
    psil=2/kxm*(sin(2/kxm));
  }

  return psil;
  }
  /* Set of Normalized thickness Distribution functions H(x) */
  double Hx_parabolic(double X)
  {
  /*
  Parabolic thickness distribution
  */
  double Hx=1-pow((2*X),2);
  return Hx;
  }

  double Hx_NACA00XX(double X)
  {
  /*
  NACA - 00XX thickness distribution
  */
  double HxNACA=5*2*(0.2969*sqrt((X+0.5))-0.1260*(X+0.5)-0.3516 \
               *pow((X+0.5),2)+0.2843*pow((X+0.5),3)-0.1015*pow((X+0.5),4));
  if(DEBUG==0) printf("x = %lf, HxNACA00XX = %lf\n",X,HxNACA);
  return HxNACA;
  }

  double Hx_NACA16seri(double X)
  {
  /*
  NACA 16 series thickness distribution - lindsey1948 (pg.5)
  */
  double HxNACA16;
  double y1=2*(0.989665*pow((X+0.5),(1/2)) - 0.239250*(X+0.5) \
               -0.041000*pow((X+0.5),2)-0.559400*pow((X+0.5),3));
  double y2=2*(0.010000+2.325000*(1-(X+0.5))-3.420000 \
               *pow((1-(X+0.5)),2) + 1.460000*pow((1-(X+0.5)),3));
  HxNACA16=MIN(y1,y2);
  return HxNACA16;
  }

  double fHx(double X)
  {
  /* choose a function to simulate */
  double Hx=Hx_NACA00XX(X);
  return Hx;
  }

  /*------------------------------------------------------------------*/
  /*          Set of Lift distribution function - f_L(x) */            
  /*------------------------------------------------------------------*/
  double fLx_uniform(double X)
  {
  /* uniform Lift distribution */
  double fL=1+0*X;
  return fL;
  }

  double fLx_parabolic(double X)
  {
  /* parabolic Lift distribution */
  double fL=3/2-6*pow(X,2);
  return fL;
  }

  double fLx(double X)
  {
  /* Lift distrubution function that chooses a function to simulate */
  double fL=fLx_parabolic(X);
  return fL;
  }

  /*------------------------------------------------------------------*/
  /*           set of Drag distribution functions                     */
  
  double fDx_uniform(double X)
  {
  /* unifrom Drag distribution */
  double fD=1+0*X;
  return fD;
  }

  double fDx_parabolic(double X)
  {
  /* parabolic Lift distribution */
  double fD=3/2-6*pow(X,2);
  return fD;
  }

  double fDx(double X)
  {
  /* Lift distrubution function that chooses a function to simulate */
  double fD=fDx_uniform(X);
  return fD;
  }

  double complex cfexp(double x,double kxx)
  {
  /*
  */
  double complex ii=(0.0+1.0*I);
  double complex expfc=(0.0+0.0*I);
  expfc=cexp(ii*kxx*x);
  return expfc;
  }

  double IntPsiV_parabolic(double x)
  {
  /* 
  */
  double fx=Hx_parabolic(x); //double complex fe=cfexp(x,kxx);
  if(DEBUG==0)printf("fx= %lf\n",fx);
  return fx;
  }

  double IntPsiV_NACA00XX(double x)
  {
  /*
  */
  double fx=Hx_NACA00XX(x);// double fe=cfexp(x,kxx);
  if(DEBUG==0)printf("x = %lf fx = %lf\n",x, fx);
  return fx;
  }

  double IntPsiL_uniform(double x)
  {
  /*
  */
  double fx=fLx_uniform(x); //double complex fe=cfexp(x,kxx);
  return fx;
  }

  double IntPsiL_parabolic(double x)
  {
  double fx=fLx_parabolic(x); //double complex fe=cfexp(x,kxx);
  return fx;
  }

  double IntPsiD_uniform(double x)
  {
  double fx=fDx_uniform(x); //double complex fe=cfexp(x,kxx);
  return fx;
  }

  double IntPsiD_parabolic(double x)
  {
  double fx=fDx_parabolic(x); //double complex fe=cfexp(x,kxx);
  return fx;
  }

  /*allocate srce plot */
  srcType allocSourcePlot(srcType src)
  {
  /*
  */
  src.kx_plt=NULL; src.kx_plt=linspace(0,20,src.plt_size);
  src.PsiV.parabolic_plt=NULL; src.PsiV.parabolic_plt=dcvector(0,src.plt_size);
  src.PsiV.NACA00XX_plt=NULL; src.PsiV.NACA00XX_plt=dcvector(0,src.plt_size);
  src.PsiL.uniform_plt=NULL; src.PsiL.uniform_plt=dcvector(0,src.plt_size);
  src.PsiL.parabolic_plt=NULL; src.PsiL.parabolic_plt=dcvector(0,src.plt_size);
  src.PsiD.uniform_plt=NULL; src.PsiD.uniform_plt=dcvector(0,src.plt_size);
  src.PsiD.parabolic_plt=NULL; src.PsiD.parabolic_plt=dcvector(0,src.plt_size);

  for(int i=0;i<src.plt_size;i++){
  //  for(int i=0;i<1;i++){
    src.PsiV.parabolic_plt[i]=cqgaus(IntPsiV_parabolic,-0.5,0.5,src.kx_plt[i]);
    src.PsiV.NACA00XX_plt[i]=cqgaus(IntPsiV_NACA00XX,-0.5,0.5,src.kx_plt[i]);
    src.PsiL.uniform_plt[i]=cqgaus(IntPsiL_uniform,-0.5,0.5,src.kx_plt[i]);
    src.PsiL.parabolic_plt[i]=cqgaus(IntPsiL_parabolic,-0.5,0.5,src.kx_plt[i]);
    src.PsiD.uniform_plt[i]=cqgaus(IntPsiD_uniform,-0.5,0.5,src.kx_plt[i]);
    src.PsiD.parabolic_plt[i]=cqgaus(IntPsiD_parabolic,-0.5,0.5,src.kx_plt[i]);

    if(DEBUG==0){
      //printf(" %lf %lf %lf\n",src.kx_plt[i],creal(src.PsiV.uniform_plt[i]),cimag(src.PsiV.uniform_plt[i]));
      //printf("kx= %lf real-naca00xx = %lf imag-naca00xx %lf\n",src.kx_plt[i],creal(src.PsiV.NACA00XX_plt[i]),cimag(src.PsiV.NACA00XX_plt[i]));
      //printf("kx= %lf real-naca00xx = %lf imag-naca00xx %lf\n",src.kx_plt[i],creal(src.PsiD.uniform_plt[i]),cimag(src.PsiD.uniform_plt[i]));

      //printf(" cexp(3) = %lf %lf\n",creal(cfexp(3,src.kx_plt[i])),cimag(cfexp(3,src.kx_plt[i])));
    }
  }
  return src;
  }

  void deallocSourcePlot(srcType src)
  {
  /*
  */
    free_dcvector(src.PsiV.parabolic_plt,0,src.plt_size);
    free_dcvector(src.PsiV.NACA00XX_plt,0,src.plt_size);
    free_dcvector(src.PsiL.uniform_plt,0,src.plt_size);
    free_dcvector(src.PsiL.parabolic_plt,0,src.plt_size);
    free_dcvector(src.PsiD.uniform_plt,0,src.plt_size);
    free_dcvector(src.PsiD.parabolic_plt,0,src.plt_size);
  return;
  }


/*File reader functions */
/*    */
  d3_t *freader3col(char* fn,int* n)
  /*-------------------------------------------------------------------
  Purpose:
             reads a data file with 3 columns and returns the data
             and number of rows. 
  Written by:
             Yunusi Fuerkaiti
             email: y.fuerkaiti@tudelft.nl
  date:
       12.06.2020
  -------------------------------------------------------------------*/
  {
  const char* thisroutine="d3_t *freader3col(...)";

  d3_t *data=NULL;
  d3_t data_tmp;
  int index,ii=0;
  const int ncol=3;
  for(int i=0;i<ncol;i++)data_tmp.x[i]=0.0;

  /* */
  FILE *file;
  if((file=fopen(fn,"r"))==NULL){
        APnoise_error("Error on opening input file.\n", thisroutine);
  }
  file=fopen(fn,"r");

  while(fscanf(file,"%lf %lf %lf",&data_tmp.x[0], \
        &data_tmp.x[1], &data_tmp.x[2])!=EOF){
    if(data==NULL){
      data=malloc(sizeof(d3_t));
      *data=data_tmp;
    }
    if(data!=NULL){
      ii++;
      data=realloc(data,sizeof(d3_t)*ii);
      index=ii-1;
      *(data+index)=data_tmp;
    }
  }
  *n=ii;
  fclose(file);
  return data;
  }
/* */

  d5_t *freader5col(char* fn,int* n)
  /*-------------------------------------------------------------------
  Purpose:
             reads a data file with 6 columns and returns the data
             and number of rows. 
  Written by:
             Yunusi Fuerkaiti
             email: y.fuerkaiti@tudelft.nl
  Date:
       12.06.2020
  Modification: 
                It reads and discards the first line that contains 
  names of each column.
  Date of modification:
  26.02.2021
  -------------------------------------------------------------------*/
  {
  const char* thisroutine="d5_t *freader5col(...)";

  d5_t *data=NULL;
  d5_t data_tmp;
  int index,ii=0;
  const int ncol=5;
  for(int i=0;i<ncol;i++) data_tmp.x[i]=0.0;
  /* */
  /* */
  FILE *file;
  if((file=fopen(fn,"r"))==NULL){
        APnoise_error("Error on opening input file.\n", thisroutine);
  }
  file=fopen(fn,"r");
  /* skip the first line*/
  fscanf(file,"%*[^\n]");//read and discard the first line
  while(fscanf(file,"%lf %lf %lf %lf %lf",&data_tmp.x[0], \
        &data_tmp.x[1], &data_tmp.x[2],&data_tmp.x[3], \
        &data_tmp.x[4])!=EOF){
    if(data==NULL){
      data=malloc(sizeof(d5_t));
      *data=data_tmp;
    }
    if(data!=NULL){
      ii++;
      data=realloc(data,sizeof(d5_t)*ii);
      index=ii-1;
      *(data+index)=data_tmp;
    }
  }
  *n=ii;
  fclose(file);
  return data;
  }


  d6_t *freader6col(char* fn,int* n)
  /*-------------------------------------------------------------------
  Purpose:
             reads a data file with 6 columns and returns the data
             and number of rows. 
  Written by:
             Yunusi Fuerkaiti
             email: y.fuerkaiti@tudelft.nl
  Date:
       12.06.2020
  Modification: 
                It reads and discards the first line that contains 
  names of each column.
  Date of modification:
  26.02.2021
  -------------------------------------------------------------------*/
  {
  const char* thisroutine="d6_t *freader6col(...)";
  d6_t *data=NULL;
  d6_t data_tmp;
  int index,ii=0;
  const int ncol=6;
  for(int i=0;i<ncol;i++) data_tmp.x[i]=0.0;
  /* */
  FILE *file;
  if((file=fopen(fn,"r"))==NULL){
        APnoise_error("Error on opening input file.\n", thisroutine);
  }
  file=fopen(fn,"r");
  fscanf(file,"%*[^\n]");//read and discard the first line
  while(fscanf(file,"%lf %lf %lf %lf %lf %lf",&data_tmp.x[0], \
        &data_tmp.x[1], &data_tmp.x[2],&data_tmp.x[3], \
        &data_tmp.x[4],&data_tmp.x[5])!=EOF){
    if(data==NULL){
      data=malloc(sizeof(d6_t));
      *data=data_tmp;
    }
    if(data!=NULL){
      ii++;
      data=realloc(data,sizeof(d6_t)*ii);
      index=ii-1;
      *(data+index)=data_tmp;
    }
  }
  *n=ii;
  fclose(file);
  return data;
  }


  d4_t *polarDataReader(char* fn,int* n)
  /*-------------------------------------------------------------------
  Purpose:
             reads a polar data file with 4 columns and returns the data
             and number of rows. 
  Written by:
             Yunusi Fuerkaiti
             email: y.fuerkaiti@tudelft.nl
  date:
       21.02.2021
  -------------------------------------------------------------------*/
  {
  const char* thisroutine="d4_t *polarDataReader(char* fn,int* n)";

  d4_t *data=NULL;
  d4_t data_tmp;
  int index,ii=0;
  const int ncol=4;
  for(int i=0;i<ncol;i++)data_tmp.x[i]=0.0;

  /* */
  FILE *file;
  if((file=fopen(fn,"r"))==NULL){
        APnoise_error("Error on opening input file.\n", thisroutine);
  }
  file=fopen(fn,"r");
  /* skipping the first two columns that start with "#" */
  fscanf(file,"%*[^\n]");
  fscanf(file,"%*[^\n]");
  while(fscanf(file,"%lf %lf %lf %lf",&data_tmp.x[0], \
        &data_tmp.x[1], &data_tmp.x[2],&data_tmp.x[3])!=EOF){
    if(data==NULL){
      data=malloc(sizeof(d4_t));
      *data=data_tmp;
    }
    if(data!=NULL){
      ii++;
      data=realloc(data,sizeof(d4_t)*ii);
      index=ii-1;
      *(data+index)=data_tmp;
    }
  }
  *n=ii;
  fclose(file);
  return data;
  }
/* */
  /*

  */
  Airfoil_t *AirfoilDBreader(char* fn,int* n)
  /*-------------------------------------------------------------------
  Purpose:
             reads airfoils data long a blade section. The first 4 
             columns are: STATION/R  CHORD/R  TWIST (DEG)  SWEEP/R 
             and the last column is airfoil PROFILE file name.
  Written by:
             Yunusi Fuerkaiti
             email: y.fuerkaiti@tudelft.nl
  date:
       19.12.2020
  -------------------------------------------------------------------*/
  {
  const char* thisroutine="Airfoil_t *AirFoilDBreader(...)";

  Airfoil_t *data=NULL;
  Airfoil_t data_tmp;
  int index,ii=0;
  const int ncol=5;
  for(int i=0;i<ncol-1;i++){
    data_tmp.x[i]=0.0;
    // data_tmp.FoilFN=" ";
  }
  /* */
  /* */
  FILE *file;
  if((file=fopen(fn,"r"))==NULL){
        APnoise_error("Error on opening input file.\n", thisroutine);
  }
  file=fopen(fn,"r");

  while(fscanf(file,"%lf %lf %lf %lf %s",&data_tmp.x[0], \
        &data_tmp.x[1],&data_tmp.x[2],&data_tmp.x[3], \
        data_tmp.FoilFN)!=EOF){
    if(data==NULL){
      data=malloc(sizeof(Airfoil_t));
      *data=data_tmp;
    }
    if(data!=NULL){
      ii++;
      data=realloc(data,sizeof(Airfoil_t)*ii);
      index=ii-1;
      *(data+index)=data_tmp;
    }
  }
  *n=ii;
  fclose(file);
  return data;
  }

/*--------------------------------------------------------------------*/
/*                                                                    */
/*                         Math utility functions                     */
/*                                                                    */
/*--------------------------------------------------------------------*/
/*--------------- Barycentric interpolation functions ----------------*/
  d2_t blii1d(d2_t x,\
              d2_t y,
              double xi)
/***********************************************************************
! Barycentric Linear Interpolation 1D
! Written by YF
! returns the interpolated value and its first derivatives:yi and yxi;
******************************************************************/
{
  double x1,x2;
  d2_t yi;
  double y1,y2;

  x1=x.x[0];
  x2=x.x[1];
  
  y1=y.x[0];
  y2=y.x[1];

  yi.x[1]=(y2-y1)/(x2-x1);
  yi.x[0]=y1+yi.x[1]*(xi-x1);
  
  return yi;
}

  int brcket(int n,\
             double x[],\
             double xi)
/* **********************************************************************
c      BRaCKET subroutine 
c      Written by YF
********************************************************************** */
{
  int i,ia,im,ib;
  
  ia = 0;
  ib = n-2;
  if(DEBUG==0) printf("%lf %lf\n",xi,x[n-1]);

/*If xi is outside the interval [x(1) x(n)] just output zero and
      return to the calling program:*/
  if((xi<x[0])||(xi>x[n-1])){
    i=0;
    if(DEBUG==0) printf("outside--brkt");
  return i;
  }
/*If xi is inside the interval [x(1) x(n)] just divide it by half
 until ib - ia = 1:*/
  im=(ia+ib)/2;
  if((ib-ia)>1){
    if((xi>x[ia]) && (xi<x[im])){
      ib=im;}
    else{
      ia=im;
    }
    im=(ia+ib)/2;
  }
  i=im;
  return i;
}

    d3_t bcui1d(d4_t x,\
                d4_t f,\
                double xi)
/*****************************************************************
! Barycentric Cubic Interpolation 1D 
! Written by Yunusi Fuerkaiti
! returns the interpolated value for given xi and its first and
! second derivates to x[0],x[1],x[2], respectively.
! 
!*****************************************************************/
{
  int i=0;
  d3_t a,px,sx,qx;
  d3_t fi;
  double x1,x2,x3,x4;
  x1=x.x[0];
  x2=x.x[1];
  x3=x.x[2];
  x4=x.x[3];

  px.x[0]= ( x2 - x1 )*( x2 - x3 )*( x2 - x4 );
  px.x[1]= ( x3 - x1 )*( x3 - x2 )*( x3 - x4 );
  px.x[2]= ( x4 - x1 )*( x4 - x2 )*( x4 - x3 );
  
  for(i=0;i<3;i++) a.x[i]=(f.x[i+1]-f.x[0])/px.x[i];
  px.x[0]=(xi-x1)*( xi - x3 )*( xi - x4 );
  px.x[1]=(xi-x1)*( xi - x2 )*( xi - x4 );
  px.x[2]=(xi-x1)*( xi - x2 )*( xi - x3 );
  
  qx.x[0]=(xi-x1)*(xi-x3)+(xi-x1)*(xi-x4)+(xi-x3)*(xi-x4);
  qx.x[1]=(xi-x1)*(xi-x2)+(xi-x1)*(xi-x4)+(xi-x2)*(xi-x4);
  qx.x[2]=(xi-x1)*(xi-x2)+(xi-x1)*(xi-x3)+(xi-x2)*(xi-x3);
  
  sx.x[0]=2.0*( 3.0*xi - x1 - x3 - x4 );
  sx.x[1]=2.0*( 3.0*xi - x1 - x2 - x4 );
  sx.x[2]=2.0*( 3.0*xi - x1 - x2 - x3 );
  
  fi.x[0]=f.x[0];//interpolated value
  fi.x[1]=0.0;//first derivative
  fi.x[2]=0.0;//second derivative
  
  for(i=0;i<3;i++){
    fi.x[0]=fi.x[0]+a.x[i]*px.x[i];
    fi.x[1]=fi.x[1]+a.x[i]*qx.x[i];
    fi.x[2]=fi.x[2]+a.x[i]*sx.x[i];}
  return fi;
}

  d3_t interp1D(int nx,\
                double tabx[],\
                double tabc[],\
                double xi)
/********************************************************************
!   (1D) Interpolation of a variable its derivatives
!   Written by YF
!********************************************************************/
{
  d3_t ci;//ci,cxi,cxxi
  d2_t c2;
  d4_t xc,fc;
  d2_t xl,fl;
  double cxxxi=0;
  int i;
  ci.x[2]=0.0;
  if(xi<tabx[0]){
    xl.x[0]=tabx[0];
    xl.x[1]=tabx[1];
    fl.x[0]=tabc[0];
    fl.x[1]=tabc[1];
    c2=blii1d(xl,fl,xi);
    ci.x[0]=c2.x[0];
    ci.x[1]=c2.x[1];
    if(DEBUG==0) printf("test1 %lf %lf\n",xi,tabx[0]);
  }
  else if(xi>tabx[nx-2]){
    xl.x[0]=tabx[nx-2];
    xl.x[1]=tabx[nx];
    fl.x[0]=tabc[nx-2];
    fl.x[1]=tabc[nx];
    c2=blii1d(xl,fl,xi);
    ci.x[0]=c2.x[0];
    ci.x[1]=c2.x[1];
    if(DEBUG==0) printf("test2 %lf %lf \n",xi,tabx[nx-2]);

  }
  else{
    i=brcket(nx,tabx,xi);
    xc.x[0]=tabx[i-1];
    xc.x[1]=tabx[i];
    xc.x[2]=tabx[i+1];
    xc.x[3]=tabx[i+2];
    
    fc.x[0]=tabc[i-1];
    fc.x[1]=tabc[i];
    fc.x[2]=tabc[i+1];
    fc.x[3]=tabc[i+2];
    ci=bcui1d(xc,fc,xi);
    if(DEBUG==0) printf("test3\n");

    }
    return ci;
}

  void ratint(double xa[],\
              double ya[],\
              int n,\
              double x,\
              double *y)
{
  char *thisroutine="void ratint(...)";
  double *dy=NULL;
  double TINY=1.0e-25;
  int m,i,ns=0;
  double w,t,hh,h,dd,*c=NULL,*d=NULL;

  c=dvector(0,n);
  d=dvector(0,n);
  hh=fabs(x-xa[0]);
  for (i=0;i<n;i++) {
    h=fabs(x-xa[i]);
    if (h == 0.0) {
            *y=ya[i];
      if(DEBUG==0) printf("y=%lf\n",*y);
             // *dy=0.0;
            //free_dvector(d,0,n);free_dvector(c,0,n);
    } else if (h < hh) {
            ns=i;
            hh=h;
    }
    c[i]=ya[i];
    d[i]=ya[i]+TINY;
  }
  *y=ya[ns--];
  for (m=1;m<n;m++) {
          for (i=0;i<=n-m;i++) {
                  w=c[i+1]-d[i];
                  h=xa[i+m]-x;
                  t=(xa[i]-x)*d[i]/h;
                  dd=t-c[i+1];
                  if (dd == 0.0) dd=dd+TINY;// APnoise_warning("Error in routine %s",thisroutine);
                  dd=w/dd;
                  d[i]=c[i+1]*dd;
                  c[i]=t*dd;
          }
          *y += (2*ns < (n-m) ? c[ns+1] : d[ns--]);
  if(DEBUG==0) printf("y=%lf\n",*y);
  }
  if(x==xa[n-1]) *y=ya[n-1];
  free_dvector(d,0,n);free_dvector(c,0,n);
}

  void polint(double xa[],\
              double ya[],\
              int n,\
              double x,\
              double *y)
  {
  char *thisroutine="void polint(...)";

        double *dy=NULL;
	int i,m,ns=1;
	double den,dif,dift,ho,hp,w;
	double *c,*d;

	dif=fabs(x-xa[0]);
	c=dvector(0,n);
	d=dvector(0,n);
	for (i=0;i<n;i++){
		if ( (dift=fabs(x-xa[i])) < dif) {
			ns=i;
			dif=dift;
		}
		c[i]=ya[i];
		d[i]=ya[i];
	}
	*y=ya[ns--];
	for (m=1;m<n;m++) {
		for (i=1;i<=n-m;i++) {
			ho=xa[i]-x;
			hp=xa[i+m]-x;
			w=c[i+1]-d[i];
			//if ( (den=ho-hp) == 0.0) APnoise_warning("Error in routine %s\n",thisroutine);
			den=w/den;
			d[i]=hp*den;
			c[i]=ho*den;
		}
		// *y += (*dy=(2*ns < (n-m) ? c[ns+1] : d[ns--]));
		*y += (2*ns < (n-m) ? c[ns+1] : d[ns--]);

	}
	free_dvector(d,0,n);
	free_dvector(c,0,n);
}


  double *diff(double *vec,\
               int N)
  /* */
  {
  double *diffvec=NULL; diffvec=dvector(0,N);
  int i;
  for(i=0;i<N;i++){
    diffvec[i]=vec[i+1]-vec[i];
  }
  return diffvec;
  }

  double *linspace(double start,\
                   double end,\
                   int N)
  /* */
  {
  double *vec=NULL;vec=dvector(0,N);
  for(int i=0;i<N;i++) vec[i]=start+(end-start)/(N-1)*i;
  return vec;
  }

  double getMax(double arr[], int n)
  {
  /* Find the maximum value in an array */
  int i;
  // Initialize maximum element
  double max = arr[0];
  // Traverse array elements from second and
  // compare every element with current max  
  for (i = 1; i < n; i++)
      if (arr[i] > max)
          max = arr[i];
  return max;
  }

  minValLoc findMinValLoc(double *arr,\
                          int numArr)
  /* */
  {
  minValLoc mVL;
  int c, location=0;

  for (c=0;c<numArr;c++)
    if (arr[c]<arr[location])
      location=c;
  mVL.minValIdx=location;
  mVL.minVal=arr[location];
  return mVL;
  }

char* cutoffstr(const char* str,\
                int from,\
                int to)
/* */
{
    if (from >= to)
        return  NULL;

    char* cut = calloc(sizeof(char), (to - from) + 1);
    char* begin = cut;
    if (!cut)
        return  NULL;
    const char* fromit = str+from;
    const char* toit = str+to;
    (void)toit;
    memcpy(cut, fromit, to);
    return begin;
}

  int strFind(Airfoil_t *profile,\
              int  N,\
              int  idx,\
              char pattern[128])
  /*
  strcmp function is a build in function that returns 0 if both strings are identical,
  for more information, please google it "strcmp()")
  */
  {
  const char* thisroutine="int strFind(...)";
  int strIdx=-1;
  char* first4chars;
  //for(int k=0;k<N;k++){
    first4chars=cutoffstr(profile[idx].FoilFN,0,4);
    if(strcmp(first4chars,pattern)==0){
      strIdx=idx;}
    if(strcmp(first4chars,pattern)!=0){
      //APnoise_warning("String can not be found.\n", thisroutine);
    }
  //}
  //printf("%s \n",first4chars);
  //printf("true %s %s\n",first4chars,pattern);
  return strIdx;
  }
/*---------------------------------------------------------------------
                       INTEGRATION FUNCTIONS                           
---------------------------------------------------------------------*/
double complex qsimp(double (*func)(double), double a, double b, double kxx)
{
  char* thisroutine="double qsimp(...)";
  double EPS=1e-6;
  int JMAX=45;
  double complex trapzd(double (*func)(double), double a, double b,double kxx, int n);
  void nrerror(char error_text[]);
  int j;
  double complex s,st,ost=(0.0+0.0*I),os=(0.0+0.0*I);

  for (j=0;j<JMAX;j++) {
          st=trapzd(func,a,b,kxx,j);
          s=(4.0*st-ost)/3.0;
          if (j > 5)
                  if (cabs(s-os) < EPS*cabs(os) ||
                          (s == 0.0 && os == 0.0)) return s;
          os=s;
          ost=st;
  }
  APnoise_error("Too many steps in routine %s",thisroutine);
  return 0.0;
}


double complex trapzd(double (*func)(double), double a, double b, double kxx, int n)
{
  char* thisroutine="double trapzd(...)";
  double complex fe1=(0.0+0.0*I);
  double x,tnm,sum,del;
  static double complex s;
  int it,j;

  if (n == 1) {
          return (s=0.5*(b-a)*(func(a)+func(b)));
  } else {
          for (it=1,j=0;j<n-1;j++) it <<= 1;
          tnm=it;
          del=(b-a)/tnm;
          x=a+0.5*del;
          for (sum=0.0,j=0;j<it;j++,x+=del){
            fe1=cfexp(x,kxx);
            sum += func(x*fe1);}
          s=0.5*(s+(b-a)*sum/tnm);
          return s;
  }
}

  double complex cqgaus(double (*func)(double), double a, double b,double kxx)
  {
    int j;
    double complex fe1=(0.0+0.0*I),fe2=(0.0+0.0*I);
    double xr,xm,dx;
    double complex s;
    static double x[]={0.0,0.1488743389,0.4333953941,
            0.6794095682,0.8650633666,0.9739065285};
    static double w[]={0.0,0.2955242247,0.2692667193,
            0.2190863625,0.1494513491,0.0666713443};

    xm=0.5*(b+a);
    xr=0.5*(b-a);
    s=(0.0+0.0*I);
    for (j=1;j<=5;j++) {
            dx=xr*x[j];
            fe1=cfexp((xm+dx),kxx);
            fe2=cfexp((xm-dx),kxx);
            s += w[j]*((*func)((xm+dx)*fe1)+(*func)((xm-dx)*fe2));
            if(DEBUG==0){
              //printf("xm+dx = %lf, xm-dx = %lf\n",xm+dx,xm-dx);
              //printf(" s= %lf %lf\n",creal(s),cimag(s));
              printf("fe1 = %lf %lf, fe2 = %lf %lf\n",creal(fe1),cimag(fe1),creal(fe2),cimag(fe2));
            }
    }
    return s *= xr;

  }


  double qgaus(double (*func)(double), double a, double b)
  {
    int j;
    double xr,xm,dx,s;
    static double x[]={0.0,0.1488743389,0.4333953941,
            0.6794095682,0.8650633666,0.9739065285};
    static double w[]={0.0,0.2955242247,0.2692667193,
            0.2190863625,0.1494513491,0.0666713443};

    xm=0.5*(b+a);
    xr=0.5*(b-a);
    s=0;
    for (j=1;j<=5;j++) {
            dx=xr*x[j];
            s += w[j]*((*func)(xm+dx)+(*func)(xm-dx));
    }
    return s *= xr;
  }

  double integTrapz(double *x, int N)
  {
  /*Purpose:
            Performs integration using the trapizoidal method
  written by:
            YF
  Date:
       23.02.2021
  */
  double integrated,sum=0.0;
  
  for(int i=0;i<N;i++){
    sum+=x[i];
  }
  integrated=(2*sum-x[0]-x[N-1])*(2*PI/N)/2;
  return integrated;
  }

  double complex integTrapzc(double complex *x, int N)
  {
  /*Purpose:
            Performs integration for a complex array using the trapizoidal method
  written by:
            YF
  Date:
       23.02.2021
  */
  double complex integrated=(0.0+0.0*I),sum=(0.0+0.0*I);
  
  for(int i=0;i<N;i++){
    sum+=x[i];
  }
  integrated=(2*sum-x[0]-x[N-1])*(2*PI/N)/2;
  return integrated;
  }

  double factorial(int n)
  {
  /* */
    int i;
    double fact=1;
    for(i=n;i>=1;i--){
      fact=fact*i;
    }
    return fact;
  }

  double besselj(int n,double x)
  {
  /* Purpose:
            Computes the Bessel function of the first kind
  Written by:
            YF
  Date:
  23.02.2021
  */
    double t0,t1,R,sum,eps;
    eps=1e-15;
        int k=1;
        //Initialize First Term
        t0=1/factorial(n);
        //Make sum equal to the first term
        sum=t0;
        do{
            //Find the ratio of the second term to the first term using already known relation
            R=-(x*x/4)/(k*(n+k));
            //Calculate the second term
            t1=R*t0;
            //find the new sum
            sum=sum+t1;
            t0=t1;
            k++;
            //keep on summing terms until the required accuracy is reached
        }while(fabs(t1/sum)>eps);
        sum=sum*pow(x/2,n);
  return sum;
}
/* Hanson's formulation functions */
  double kx(PropellerData prop,\
            double c_R,\
            double r_R,\
            double theta_r,\
            int m)
  {
  /*
  */
  double kx_,Mach_r;
  Mach_r=sqrt(prop.Mach_x*prop.Mach_x+r_R*prop.Mach_t*prop.Mach_t);
  kx_=(2*m*prop.NB*c_R)/(Mach_r*(1-prop.Mach_x*cos(theta_r)));
  
  return kx_;
  }

  double Psi_L(double kx)
  {
  /* 
  */
  double psiL;
  if(kx==0){
    psiL=1;}
  else{
    psiL=2/kx*(sin(kx/2));
  }
  return psiL;
  }


  double get_cpu_time(clock_t start,\
                      clock_t end)
  {
  /*
  */
  double cpu_time_used=0.0;
  cpu_time_used=((double) (end-start))/CLOCKS_PER_SEC;
  printf("\n\n");
  printf(" -----------------------------------------------\n");
  printf(" Computation completed sucessfully! \n");
  printf(" CPU time %lf seconds\n", cpu_time_used);
  printf(" -----------------------------------------------\n");
  return cpu_time_used;
}


