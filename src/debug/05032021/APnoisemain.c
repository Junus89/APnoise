#include "APnoise.h"



int main(int argc, char **argv)
{
  /*cpu time variables */
  clock_t start=0.0;
  clock_t end=0.0;
  double cpu_time;
  start=clock();
  /*reading input file */
  const char* thisroutine="int main(...)";
  /* input data declerations */
  caseData        caseIN;
  PropellerData   prop;
  AtmData         atm;
  ObserverData    obsrvr;
  char skipLine[256];
  FILE *fp=NULL;
  if(argc>=2){
    fp=fopen(argv[1],"r");
  }
  else{
    APnoise_error("Can no read the input file.\n\
Please check your input file!\n",thisroutine);
  }
  while(1==fscanf(fp,"%s%*[^\n] \
                      %d%*[^\n] \
                      %s%*[^\n] \
                      %lf%*[^\n] \
                      %lf%*[^\n] \
                      %lf%*[^\n] \
                      %d%*[^\n] \
                      %lf%*[^\n] \
                      %lf%*[^\n] \
                      %lf%*[^\n] \
                      %lf%*[^\n] \
                      %s%*[^\n] \
                      %s%*[^\n] \
                      %s%*[^\n] \
                      %lf%*[^\n] \
                      %lf%*[^\n] \
                      %lf%*[^\n] \
                      %lf%*[^\n] \
                      %lf%*[^\n] \
                      %s%*[^\n] \
                      %d%*[^\n] \
                      %s%*[^\n] \
                      %d%*[^\n] \
                      %d%*[^\n] \
                      %s%*[^\n] \
                      %s%*[^\n] \n",\
                      caseIN.caseName,\
                      &caseIN.runID,\
                      skipLine,\
                      &prop.D,\
                      &prop.bladeAngle,\
                      &prop.bladeChord,\
                      &prop.NB,\
                      &prop.RPM,\
                      &prop.Mach_x,\
                      &prop.alpha,\
                      &prop.CP,\
                      prop.bladeGeom.Fn,\
                      prop.bladeAero.Fn,\
                      skipLine,\
                      &atm.alt,\
                      &atm.T,\
                      &atm.p,\
                      &atm.rho,\
                      &atm.c0,\
                      skipLine,\
                      &obsrvr.ifRead,\
                      obsrvr.coordsFN,\
                      &obsrvr.thetaNum,\
                      &obsrvr.phiNum,\
                      skipLine,\
                      prop.DbDev.FN)){};
/*check if we get the correct input data */
prop.J=(prop.Mach_x*atm.c0)/(prop.RPM/60*prop.D);
prop.eta=prop.J*fabs((prop.CT/prop.CP));
prop.omega=prop.RPM/30.0*PI;
prop.Mach_t=prop.omega*prop.D/2.0/atm.c0;
prop.Mach_h=sqrt(pow(prop.Mach_t,2)+pow(prop.Mach_x,2));
 prop.Re = 0.95*prop.D/2;  /* effective radius */
prop.HNum=10;


double p_ref=2e-5;
if(DEBUG==1){
  printf("============================================\n\n");
  printf("           DEBUG INFO ENABLED               \n\n");
  printf("case name %s and NB %d and CT  %lf\n",caseIN.caseName,prop.NB,prop.CT);
  printf("Blade geom data fn: %s\n",prop.bladeGeom.Fn);
  printf("Blade aero data fn: %s\n",prop.bladeAero.Fn);
  printf(" reference pressure %lf\n",p_ref);
  printf("ThetaNum = %d, and PhiNum = %d\n",obsrvr.thetaNum,obsrvr.phiNum);
  printf("Omega = %lf\n",prop.omega);
  printf("Eta =%lf\n",prop.eta);
  printf("Mach_h = %lf, Mach_tip = %lf\n",prop.Mach_h,prop.Mach_t);
  printf("J = %lf\n",prop.J);
}

/* read or generate receivers */
int FFobNum=0;
int XNum=obsrvr.thetaNum;
int ZNum=1;
int FFobsNumAll=XNum*ZNum;
//double *FFx=NULL;
double *FFx=NULL;
double *FFxAngle=NULL;
double *FFy=NULL;
double *FFz=NULL;


if(obsrvr.ifRead==1){
  /* read obsrver geometry */
  obsrvr.FFcoords=freader3col(obsrvr.coordsFN,&obsrvr.FFNum);
  obsrvr.Y=NULL;obsrvr.Y=dvector(0,obsrvr.FFNum);
  obsrvr.S0=NULL;obsrvr.S0=dvector(0,obsrvr.FFNum);
  obsrvr.S=NULL;obsrvr.S=dvector(0,obsrvr.FFNum);
  obsrvr.theta=NULL;obsrvr.theta=dvector(0,obsrvr.FFNum);

  double Bet=sqrt(1-prop.Mach_x*prop.Mach_x);

  for(int i=0;i<obsrvr.FFNum;i++){
    obsrvr.S[i]=sqrt(pow(obsrvr.FFcoords[i].x[0],2)+pow(obsrvr.FFcoords[i].x[1],2)+pow(obsrvr.FFcoords[i].x[2],2));
    obsrvr.Y[i]=sqrt(pow(obsrvr.FFcoords[i].x[1],2)+pow(obsrvr.FFcoords[i].x[2],2));
    obsrvr.S0[i]=sqrt(pow(obsrvr.FFcoords[i].x[0],2)+Bet*Bet*pow(obsrvr.Y[i],2));
    obsrvr.theta[i]=acos(obsrvr.FFcoords[i].x[0]/obsrvr.S[i]);
  }
  FFxAngle=dvector(0,obsrvr.FFNum);
  for(int i=0;i<obsrvr.FFNum;i++){
    FFxAngle[i]=180-atan2(obsrvr.FFcoords[i].x[1],obsrvr.FFcoords[i].x[0])*180/PI;
    if(DEBUG==1){
      printf(" x = %lf, y = %lf, theta = %lf\n",obsrvr.FFcoords[i].x[0],\
                      obsrvr.FFcoords[i].x[1],FFxAngle[i]);
      printf(" obsrvr.theta[%d] = %lf\n",i,(180-obsrvr.theta[i]*180/PI));
    }
  }
  if(DEBUG==1){
   printf("Receiver file readed:obsrvr.FFNum %d\n",obsrvr.FFNum);
   printf("Atan2(1,1) = %lf, atan(1) = %lf\n",atan2(PI/4,1),atan(PI/4));
  }
}
else{
FFobNum=0;
XNum=obsrvr.thetaNum;
ZNum=1;
FFobsNumAll=XNum*ZNum;
FFx=linspace(0,PI,XNum);FFxAngle=linspace(0,PI,XNum);
for(int i=0;i<XNum;i++) FFx[i]=0.7*prop.D*cos(FFx[i]);

FFy=linspace(-50*prop.D,50*prop.D,XNum);
FFz=linspace(-10*prop.D,10*prop.D,ZNum);
obsrvr.FFNum=0; obsrvr.FFcoords=NULL; obsrvr.FFcoords=d3_tvector(0,FFobsNumAll);
for(int i=0;i<XNum;i++){
  for(int j=0;j<ZNum;j++){
    obsrvr.FFNum++;
    obsrvr.FFcoords[obsrvr.FFNum-1].x[0]=FFx[obsrvr.FFNum-1];
    obsrvr.FFcoords[obsrvr.FFNum-1].x[1]=0.8*prop.D;
    obsrvr.FFcoords[obsrvr.FFNum-1].x[2]=0.0;
    if(DEBUG==0) printf("idx %d, x = %lf\n",obsrvr.FFNum,obsrvr.FFcoords[obsrvr.FFNum-1].x[0]);
  }
}
}
if(DEBUG==1)printf("Number of farfield obsrvers %d\n",obsrvr.FFNum);
/* 
*/
//prop.Thrust = prop.CT*atm.rho*(pow((double)prop.RPM/60,2)*pow(prop.D,4)); /* cruise thrust */
prop.Thrust = prop.CT*atm.rho*(pow((double)prop.RPM/60.0,2.0)*pow(prop.D,4));
prop.Torque = prop.Thrust*prop.Mach_x*atm.c0/prop.eta/prop.omega; /*cruise torque */
int AZNum = 3601;    /* azimuthal number for blade */
if(DEBUG==1)printf("Propeller thrust %lf\n",prop.Thrust);
if(DEBUG==1)printf("Propeller torque %lf\n",prop.Torque);
if(DEBUG==1){
  printf("\n\n");
  printf(" 6Mic coords: %lf %lf %lf \n",obsrvr.FFcoords[0].x[0],obsrvr.FFcoords[0].x[1],obsrvr.FFcoords[0].x[2]);
}

/* 
*/
double *srcTheta=NULL; srcTheta=linspace(0,TWOPI,AZNum);
d3_t *srcPos=NULL; srcPos=d3_tvector(0,AZNum);
double Beta,*S=NULL,*Sigma=NULL,kc,*AS=NULL,*BS=NULL;
S=dvector(0,AZNum);Sigma=dvector(0,AZNum);
AS=dvector(0,AZNum);BS=dvector(0,AZNum);
double A=0,B=0;
double **prms1=NULL; prms1=dmatrix(0,prop.HNum,0,obsrvr.FFNum);
double **SPL1=NULL; SPL1=dmatrix(0,prop.HNum,0,obsrvr.FFNum);
enum runtype run_t;
run_t=caseIN.runID;
FILE *fspl1;
double S0,freq,Y;
freq=prop.NB*prop.omega/2/PI;
double **prms2=NULL; prms2=dmatrix(0,prop.HNum,0,obsrvr.FFNum);
double **SPL2=NULL; SPL2=dmatrix(0,prop.HNum,0,obsrvr.FFNum);
FILE *fspl2;
/* */
FILE *fspl,*fsplm,*fsplr,*fsplr2m;
double **prmsL=NULL; prmsL=dmatrix(0,prop.HNum,0,obsrvr.FFNum);
double **prmsT=NULL; prmsT=dmatrix(0,prop.HNum,0,obsrvr.FFNum);
double **SPL=NULL; SPL=dmatrix(0,prop.HNum,0,obsrvr.FFNum);
double **SPLL=NULL; SPLL=dmatrix(0,prop.HNum,0,obsrvr.FFNum);
double **SPLT=NULL; SPLT=dmatrix(0,prop.HNum,0,obsrvr.FFNum);

switch(run_t){
  case(GutinDeming):
    fspl2=fopen("SPLGD.txt","w");
    char fnc[256];strcpy(fnc,caseIN.caseName);
    strcat(fnc,"_SPLGDm.txt");
    fsplm=fopen(fnc,"w");
    char fnr1[256];strcpy(fnr1,caseIN.caseName);
    strcat(fnr1,"_SPLGDr.txt");fsplr=fopen(fnr1,"w");
    char fnr1k[256];strcpy(fnr1k,caseIN.caseName);
    strcat(fnr1k,"_SPLGDr2m.txt");fsplr2m=fopen(fnr1k,"w");
    for(int m=1;m<prop.HNum+1;m++){
      for(int j=0;j<obsrvr.FFNum;j++){
        Beta=sqrt(1-prop.Mach_x*prop.Mach_x);
        Y=sqrt(pow(obsrvr.FFcoords[j].x[1],2)+pow(obsrvr.FFcoords[j].x[2],2));
        S0=sqrt(pow(obsrvr.FFcoords[j].x[0],2)+Beta*Beta*pow(Y,2));
        kc=m*prop.NB*prop.omega/atm.c0;

        prms2[m][j]=m*prop.NB*prop.omega/(2*PI*sqrt(2)*atm.c0*S0) \
                    *fabs(prop.Thrust*(prop.Mach_x+obsrvr.FFcoords[j].x[0]/S0)/Beta*Beta \
                    -prop.Torque*prop.NB*atm.c0/(prop.NB*prop.omega*pow(prop.Re,2))) \
                    *besselj(m*prop.NB,kc*prop.Re*obsrvr.FFcoords[j].x[1]/S0);
        if(prms2[m][j]<0.0){
          prms2[m][j]=-prms2[m][j];
          SPL2[m][j]=20*log10(prms2[m][j]/p_ref);
        }
        else{
          SPL2[m][j]=20*log10(prms2[m][j]/p_ref);
        }
        fprintf(fspl2,"%lf %lf\n",m*freq,SPL2[m][j]);
//        if(m==1) fprintf(fsplr,"%lf %lf\n",(obsrvr.FFcoords[j].x[0]/(prop.D)),SPL2[1][j]);
        if(m==1) fprintf(fsplr,"%lf %lf\n",(180-obsrvr.theta[j]*180.0/PI),SPL2[1][j]);
        if(m==2) fprintf(fsplr2m,"%lf %lf\n",(180-obsrvr.theta[j]*180.0/PI),SPL2[2][j]);
      }
      fprintf(fsplm,"%lf %lf \n",m*freq,SPL2[8][0]);
    }
    fclose(fspl2);
    fclose(fsplm);
    fclose(fsplr);
    /* second formulation */
    if(DEBUG==1){
      printf("kc, %lf, S0 = %lf\n",kc,S0);
      printf("prms2[0][0] = %lf, prms2[HNum][obsFFNum] = %lf\n",prms2[1][0],prms2[prop.HNum][obsrvr.FFNum-1]);
    }
    break;
  case(GarrickWatkins):
    fspl1=fopen("SPLGW.txt","w");
    char fn[256];strcpy(fn,caseIN.caseName);
    strcat(fn,"_SPLGWm.txt");
    fsplm=fopen(fn,"w");
    char fnr2[256];strcpy(fnr2,caseIN.caseName);
    strcat(fnr2,"_SPLGWr.txt");fsplr=fopen(fnr2,"w");
    char fnr2k[256];strcpy(fnr2k,caseIN.caseName);
    strcat(fnr2k,"_SPLGWr2m.txt");fsplr2m=fopen(fnr2k,"w");
    for(int m=1;m<prop.HNum+1;m++){
      for(int j=0;j<obsrvr.FFNum;j++){
        /* source position over time */
        for(int k=0;k<AZNum;k++){
          srcPos[k].x[0]=0.0;
          srcPos[k].x[1]=prop.Re*cos(srcTheta[k]);
          srcPos[k].x[2]=prop.Re*sin(srcTheta[k]);
          Beta=sqrt((1-pow(prop.Mach_x,2)));
          S[k]=sqrt(pow((obsrvr.FFcoords[j].x[0]-srcPos[k].x[0]),2)+
               pow(Beta,2)*(pow((obsrvr.FFcoords[j].x[1]-srcPos[k].x[1]),2)+
               pow((obsrvr.FFcoords[j].x[2]-srcPos[k].x[2]),2)));
            
          Sigma[k]=(prop.Mach_x*(obsrvr.FFcoords[j].x[0]-srcPos[k].x[0])+S[k])/pow(Beta,2);
          kc=m*prop.NB*prop.omega/atm.c0;
        /* */
        AS[k]=(prop.Thrust*obsrvr.FFcoords[j].x[0]*cos((m+1)*prop.NB*srcTheta[k]+kc*Sigma[k])/pow(S[k],2) \
             +(prop.Thrust*kc*(prop.Mach_x+obsrvr.FFcoords[j].x[0]/S[k])/pow(Beta,2) \
              -prop.Torque*m*prop.NB/pow(prop.Re,2))*sin(m*prop.NB*srcTheta[k]+kc*Sigma[k]))/S[k];

        BS[k]=(-prop.Thrust*obsrvr.FFcoords[j].x[0]*sin((m+1)*prop.NB*srcTheta[k]+kc*Sigma[k])/pow(S[k],2) \
             +(prop.Thrust*kc*(prop.Mach_x+obsrvr.FFcoords[j].x[0]/S[k])/pow(Beta,2) \
              -prop.Torque*m*prop.NB/pow(prop.Re,2))*cos(m*prop.NB*srcTheta[k]+kc*Sigma[k]))/S[k];
        }
        A=integTrapz(AS,AZNum); B=integTrapz(BS,AZNum);
        prms1[m][j]=sqrt(2)*sqrt(A*A+B*B)/(8*PI*PI);
        SPL1[m][j]=20*log10(prms1[m][j]/p_ref);
        fprintf(fspl1,"%d %d %lf\n",m,j,SPL1[m][j]);
        if(m==1) fprintf(fsplr,"%lf %lf\n",(180-obsrvr.theta[j]*180.0/PI),SPL1[1][j]);
        if(m==2) fprintf(fsplr2m,"%lf %lf\n",(180-obsrvr.theta[j]*180.0/PI),SPL1[2][j]);

      }
      fprintf(fsplm,"%lf %lf \n",m*freq,SPL1[m][8]);
    }
    fclose(fspl1);
    fclose(fsplm);
    fclose(fsplr);

    if(DEBUG==1){
      printf("kc = %lf\n",kc);
      printf("Beta = %lf\n",Beta);
      printf("Omega = %lf\n",prop.omega);
      printf("Thrust = %lf, CT = %lf\n",prop.Thrust,prop.CT);
      printf("Torque = %lf\n",prop.Torque);
      printf("RPM/60 = %lf \n",(double)(prop.RPM/60));
      printf("srctheta[0] = %lf, srctheta[AZNum-1] %lf\n",srcTheta[0],srcTheta[AZNum-1]);
      printf("S[0] = %lf, S[AZNum] = %lf\n",S[0],S[AZNum-1]);
      printf("Sigma[0] = %lf, Sigma[AZNum] = %lf\n",Sigma[0],Sigma[AZNum-1]);
      printf("AS[0] = %lf, AS[AZNum] = %lf\n",AS[0],AS[AZNum-1]);
      printf("Re = %lf \n",prop.Re);
      printf("A = %lf, B = %lf\n",A,B);
      FILE *fpAS=NULL;
      fpAS=fopen("AS.txt","w");
      for(int i=0;i<AZNum;i++) fprintf(fpAS,"%d %lf\n",i,AS[i]);
      fclose(fpAS);
    }
    break;
  case(BarryMagliozzi):

    prop=allocDbDevInput(prop);
    double integRoot2TipL=0.0,integRoot2TipT=0.0;
    double *integRTL=NULL; integRTL=dvector(0,prop.DbDev.NumStations);
    double *integRTT=NULL; integRTT=dvector(0,prop.DbDev.NumStations);

    fspl=fopen("SPLBM.txt","w");
    char fnm[256];strcpy(fnm,caseIN.caseName);
    char fnr[256];strcpy(fnr,caseIN.caseName);
    char fnr2m[256];strcpy(fnr2m,caseIN.caseName);

    strcat(fnm,"_SPLBMm.txt");
    strcat(fnr,"_SPLBMr.txt");
    strcat(fnr2m,"_SPLBMr2m.txt");
    fsplm=fopen(fnm,"w");
    fsplr=fopen(fnr,"w");
    fsplr2m=fopen(fnr2m,"w");
    
    for(int m=1;m<prop.HNum+1;m++){
      for(int j=0;j<obsrvr.FFNum;j++){
        Beta = sqrt(1-prop.Mach_x*prop.Mach_x);
        Y=sqrt(pow(obsrvr.FFcoords[j].x[1],2)+pow(obsrvr.FFcoords[j].x[2],2));
        S0=sqrt(pow(obsrvr.FFcoords[j].x[0],2)+Beta*Beta*pow(Y,2));
        //S0 = sqrt(pow(obsrvr.FFcoords[j].x[0],2)+Beta*Beta*pow(obsrvr.FFcoords[j].x[1],2));
        kc=m*prop.NB*prop.omega/atm.c0;
        for(int k=0;k<prop.DbDev.NumStations;k++){
          prop.DbDev.Bangle[k]=prop.DbDev.Bangle[k]*PI/180.0;
          /* Loading term */
          integRTL[k]=prop.DbDev.r_R[k]/(prop.DbDev.Bchord[k]*cos(prop.DbDev.Bangle[k])) \
                     *sin(m*prop.NB*prop.DbDev.Bchord[k]*cos(prop.DbDev.Bangle[k])/(2*prop.DbDev.r_R[k])) \
                     *((prop.Mach_x+obsrvr.FFcoords[j].x[0]/S0)*prop.omega*prop.DbDev.dT_dr[k]/(atm.c0*Beta*Beta) \
                     -1/(prop.DbDev.r_R[k]*prop.DbDev.r_R[k])*prop.DbDev.dQ_dr[k]) \
                     *(besselj(m*prop.NB,kc*prop.DbDev.r_R[k]*obsrvr.FFcoords[j].x[1]/S0) \
                     +(pow(Beta,2)*obsrvr.FFcoords[j].x[1]*prop.DbDev.r_R[k])/(2*S0*S0) \
                     *(besselj((m*prop.NB-1),kc*prop.DbDev.r_R[k]*obsrvr.FFcoords[j].x[1]/S0) \
                     -besselj((m*prop.NB+1),kc*prop.DbDev.r_R[k]*obsrvr.FFcoords[j].x[1]/S0)));
          /* Thickness term */
          integRTT[k]=-prop.DbDev.AA[k]*(besselj(m*prop.NB,kc*prop.DbDev.r_R[k]*obsrvr.FFcoords[j].x[1]/S0) \
                      +(pow(Beta,2)*obsrvr.FFcoords[j].x[1]*prop.DbDev.r_R[k])/(2*S0*S0) \
                      *(besselj((m*prop.NB-1),kc*prop.DbDev.r_R[k]*obsrvr.FFcoords[j].x[1]/S0) \
                      -besselj((m*prop.NB+1),kc*prop.DbDev.r_R[k]*obsrvr.FFcoords[j].x[1]/S0)));

        }
        integRoot2TipL=integTrapz(integRTL,prop.DbDev.NumStations);
        prmsL[m][j]=1.0/(sqrt(2)*PI*S0)*integRoot2TipL;

        integRoot2TipT=integTrapz(integRTT,prop.DbDev.NumStations);
        prmsT[m][j]=atm.rho*pow(m,2)*pow(prop.omega,2)*pow(prop.NB,3)/(2*sqrt(2)*PI*pow(Beta,4)) \
                   *(pow((S0+prop.Mach_x*obsrvr.FFcoords[j].x[0]),2)/pow(S0,3))*integRoot2TipT;

        SPLL[m][j]=20*log10(fabs(prmsL[m][j])/p_ref);
        SPLT[m][j]=20*log10(fabs(prmsT[m][j])/p_ref);
        SPL[m][j]=20*log10(fabs(prmsL[m][j]+prmsT[m][j])/p_ref);

        fprintf(fspl,"%d %d %lf\n",m,j,SPL[m][j]);
        if(m==1) fprintf(fsplr,"%4.2f %4.4f %4.4f %4.4f\n",(180-obsrvr.theta[j]*180.0/PI),SPLT[1][j],SPLL[1][j],SPL[1][j]);
        if(m==2) fprintf(fsplr2m,"%lf %lf\n",(180-obsrvr.theta[j]*180.0/PI),SPL[2][j]);

      }
      fprintf(fsplm,"%lf %lf %lf %lf\n",m*freq,SPLT[m][0],SPLL[m][0],SPL[m][0]);
    }
    fclose(fspl);
    fclose(fsplm);
    fclose(fsplr);
    fclose(fsplr2m);
    deallocDbDevInput(&prop);
    free_dvector(integRTL,0,prop.DbDev.NumStations);
    free_dvector(integRTT,0,prop.DbDev.NumStations);
    printf("%lf\n",prop.DbDev.r_R[prop.DbDev.NumStations-10]);

    break;
  case(Hanson):

    printf("\n\nHanson's method\n\n");
    /*reading source input data */
    prop=allocDbDevInput(prop);
    /*getting retarded time varibales */
    obsrvr=getRetard(prop, obsrvr);
    /* */
    getNoise(caseIN,prop,atm,obsrvr);
/*
    prop.Mach_t=0.2308;// * try to have this from other dependent variables * /
    int N; N=100; // * number of stations of the blade * /
    d5_t *BladeGeom=NULL; BladeGeom=freader5col(prop.bladeGeom.Fn,&prop.bladeGeom.NumStations);
    prop.bladeGeom.thicknsRatio=0.12;// * needs to be read from input data? * //
    prop=allocBladeGeom(prop, BladeGeom,N);
    // * Blade Aero data * /
    d6_t *BladeAero=NULL; BladeAero=freader6col(prop.bladeAero.Fn,&prop.bladeAero.NumStations);
    prop=allocBladeAeroData(prop, BladeAero,N);
    if(DEBUG==0){
      printf("Num Stations  %d\n",prop.bladeGeom.NumStations);
      printf("Num Stations  %d\n",prop.bladeAero.NumStations);
      printf("r_R[0] = %lf\n",prop.bladeGeom.r_R[0]);
      for(int i=0;i<prop.bladeGeom.NumStations;i++){
        printf(" r_R[%d] = %lf, Cl[%d] = %lf thickness = %lf\n",i,prop.bladeGeom.r_R[i],i,prop.bladeAero.Cl[i],prop.bladeGeom.thickns[i]);
        printf(" c_D[%d] = %lf, lambda[%d] = %lf\n",i,prop.bladeGeom.c_D[i],i,prop.bladeGeom.lambdaManufac[i]);
      }
    }
    // * plot distributions * /
    prop.atZ.NumStations=N;
    plotDistributions(prop,N);
    // * source allocation * /
    srcType src; src.plt_size=201;
    src=allocSourcePlot(src);
    // * compute the noise harmonics * /
    ListenerMics mics; mics.numMics=100; mics.numHarmo=4;
    mics=getNoiseHarmonics(src,prop,atm,mics);






    free(BladeGeom);
    deallocBladeGeom(prop,N);
    //free(BladeAero);
    //deallocBladeAeroData(prop,N);
    //deallocSourcePlot(src);
    deallocNoiseMics(mics);
*/
    deallocDbDevInput(&prop);
    break;
  default:
    APnoise_error("Please specify the run type\n",thisroutine);
}
  if(DEBUG==1){
  printf("                                            \n");
  printf("           End of DEBUG INFO                \n");
  printf("============================================\n\n");
  }

free_dvector(S,0,AZNum);
free_dvector(Sigma,0,AZNum);
free_dvector(AS,0,AZNum);
free_dvector(BS,0,AZNum);
free_d3_tvector(obsrvr.FFcoords,0,FFobsNumAll);
free_d3_tvector(srcPos,0,AZNum);
free_dmatrix(prms1,0,prop.HNum,0,AZNum);
free_dmatrix(SPL1,0,prop.HNum,0,AZNum);
free_dmatrix(prms2,0,prop.HNum,0,AZNum);
free_dmatrix(SPL2,0,prop.HNum,0,AZNum);
free_dmatrix(prmsL,0,prop.HNum,0,AZNum);
free_dmatrix(prmsT,0,prop.HNum,0,AZNum);
free_dmatrix(SPL,0,prop.HNum,0,AZNum);
free_dmatrix(SPLL,0,prop.HNum,0,AZNum);
free_dmatrix(SPLT,0,prop.HNum,0,AZNum);
if(obsrvr.ifRead==1){
//  free_d3_tvector(obsrvr.FFcoords,0,obsrvr.FFNum);
  free_dvector(obsrvr.Y,0,obsrvr.FFNum);
  free_dvector(obsrvr.S,0,obsrvr.FFNum);
  free_dvector(obsrvr.S0,0,obsrvr.FFNum);
}
/* computing the CPU time */
end=clock();
cpu_time=get_cpu_time(start,end);

return 0;
}


