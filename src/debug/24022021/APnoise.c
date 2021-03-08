#include "APnoise.h"

/* utility functions */
  void APnoise_error(char error_text[],\
                     const char* routine)
  /*-------------------------------------------------------------------
  Purpose:
          reports error for the BEMT Solver implementation.

  Written by:
             Yunusi Fuerkaiti
             email: y.fuerkaiti@tudelft.nl
  date:
       12.06.2020
  -------------------------------------------------------------------*/
  {
    fprintf(stderr,"                                    \n");
    fprintf(stderr,"BEMT run-time error in routine:\t%s\n",\
      routine);
    fprintf(stderr,"%s\n",error_text);
    fprintf(stderr,"Now exiting to system ...\n");
    exit(1);
  }

  void APnoise_warning(char error_text[],\
                       const char* routine)
  /*-------------------------------------------------------------------
  Purpose:
          reports warning for the BEMT Solver implementation.

  Written by:
             Yunusi Fuerkaiti
             email: y.fuerkaiti@tudelft.nl
  date:
       12.06.2020
  -------------------------------------------------------------------*/
  {
    fprintf(stderr,"                                    \n");
    fprintf(stderr,"BEMT run-time warning in routine:\t%s\n",\
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
  printf("%lf %lf\n",xi,x[n-1]);

/*If xi is outside the interval [x(1) x(n)] just output zero and
      return to the calling program:*/
  if((xi<x[0])||(xi>x[n-1])){
    i=0;
    printf("outside--brkt");
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
!   CVALS1.FOR
!   (1D) Interpolation of sound speed and its derivatives
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
    if(DEBUG==1) printf("test1 %lf %lf\n",xi,tabx[0]);
  }
  else if(xi>tabx[nx-2]){
    xl.x[0]=tabx[nx-2];
    xl.x[1]=tabx[nx];
    fl.x[0]=tabc[nx-2];
    fl.x[1]=tabc[nx];
    c2=blii1d(xl,fl,xi);
    ci.x[0]=c2.x[0];
    ci.x[1]=c2.x[1];
    if(DEBUG==1) printf("test2 %lf %lf \n",xi,tabx[nx-2]);

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
    if(DEBUG==1) printf("test3\n");

    }
    return ci;
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
      //BEMT_warning("String can not be found.\n", thisroutine);
    }
  //}
  //printf("%s \n",first4chars);
  //printf("true %s %s\n",first4chars,pattern);
  return strIdx;
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
  integrated=(2*sum-x[0]-x[N])*(2*PI/N)/2;
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


