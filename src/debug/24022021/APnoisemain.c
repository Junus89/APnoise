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
                      %d%*[^\n] \
                      %lf%*[^\n] \
                      %lf%*[^\n] \
                      %lf%*[^\n] \
                      %s%*[^\n] \
                      %lf%*[^\n] \
                      %lf%*[^\n] \
                      %lf%*[^\n] \
                      %lf%*[^\n] \
                      %lf%*[^\n] \
                      %s%*[^\n] \
                      %d%*[^\n] \
                      %d%*[^\n] \n",\
                      caseIN.caseName,\
                      &caseIN.runID,\
                      skipLine,\
                      &prop.D,\
                      &prop.bladeAngle,\
                      &prop.bladeChord,\
                      &prop.NB,\
                      &prop.RPM,\
                      &prop.Mach,\
                      &prop.CT,\
                      &prop.eta,\
                      skipLine,\
                      &atm.alt,\
                      &atm.T,\
                      &atm.p,\
                      &atm.rho,\
                      &atm.c0,\
                      skipLine,\
                      &obsrvr.thetaNum,\
                      &obsrvr.phiNum)){};
/*check if we get the correct input data */
prop.eta=prop.eta/100.0; /* convert % to numeric */
prop.omega=prop.RPM/30.0*PI;
double p_ref=2e-5;
if(DEBUG==1){
  printf("case name %s and NB %d and CT  %lf\n",caseIN.caseName,\
  prop.NB,prop.CT);
  printf(" reference pressure %lf\n",p_ref);
  printf("ThetaNum = %d, and PhiNum = %d\n",obsrvr.thetaNum,obsrvr.phiNum);
  printf("Omega = %lf\n",prop.omega);
}
/* nearfield observers */
obsrvr.theta=NULL; obsrvr.theta=linspace(0,PI,obsrvr.thetaNum);
obsrvr.phi=NULL; obsrvr.phi=linspace(0,TWOPI,obsrvr.phiNum);
obsrvr.NFNum=0;
int obNumAll=obsrvr.thetaNum*obsrvr.phiNum;
obsrvr.NFcoords=NULL;obsrvr.NFcoords=d3_tvector(0,obNumAll);
for(int i=0;i<obsrvr.thetaNum;i++){
  for(int j=0;j<obsrvr.phiNum;j++){
    obsrvr.NFNum++;
    obsrvr.NFcoords[obsrvr.NFNum].x[0]=2*prop.D*cos(obsrvr.theta[i]);
    obsrvr.NFcoords[obsrvr.NFNum].x[1]=2*prop.D*sin(obsrvr.theta[i])*cos(obsrvr.phi[j]);
    obsrvr.NFcoords[obsrvr.NFNum].x[2]=2*prop.D*sin(obsrvr.theta[i])*sin(obsrvr.phi[j]);
  }
}
if(DEBUG==1)printf("Number of nearfield obsrvers %d\n",obsrvr.NFNum);


int FFobNum=0;
int XNum=obsrvr.thetaNum;
int ZNum=1;
int FFobsNumAll=XNum*ZNum;
//double *FFx=NULL;FFx=linspace(-50*prop.D,50*prop.D,XNum);
double *FFx=NULL;FFx=linspace(-PI*5.2,PI*5.2*prop.D,XNum);

double *FFy=NULL;FFy=linspace(-50*prop.D,50*prop.D,XNum);
double *FFz=NULL;FFz=linspace(-10*prop.D,10*prop.D,ZNum);
obsrvr.FFNum=0; obsrvr.FFcoords=NULL; obsrvr.FFcoords=d3_tvector(0,FFobsNumAll);
for(int i=0;i<XNum;i++){
  for(int j=0;j<ZNum;j++){
    obsrvr.FFNum++;
    obsrvr.FFcoords[obsrvr.FFNum-1].x[0]=FFx[obsrvr.FFNum-1];
    obsrvr.FFcoords[obsrvr.FFNum-1].x[1]=-100.0*prop.D;
    obsrvr.FFcoords[obsrvr.FFNum-1].x[2]=0.0;
    if(DEBUG==1) printf("idx %d, x = %lf\n",obsrvr.FFNum,obsrvr.FFcoords[obsrvr.FFNum-1].x[0]);

  }
}
if(DEBUG==1)printf("Number of farfield obsrvers %d\n",obsrvr.FFNum);
/* 
*/
//prop.Thrust = prop.CT*atm.rho*(pow((double)prop.RPM/60,2)*pow(prop.D,4)); /* cruise thrust */
prop.Thrust = prop.CT*atm.rho*(pow((double)prop.RPM/60.0,2.0)*pow(prop.D,4));
prop.Torque = prop.Thrust*prop.Mach*atm.c0/prop.eta/prop.omega; /*cruise torque */
prop.Re = 0.8*prop.D/2;  /* effective radius */
int HNum = 4;   /* harmonic number */
int AZNum = 3601;    /* azimuthal number for blade */
if(DEBUG==1)printf("Propeller thrust %lf\n",prop.Thrust);
/* 
*/
double *srcTheta=NULL; srcTheta=linspace(0,TWOPI,AZNum);
d3_t *srcPos=NULL; srcPos=d3_tvector(0,AZNum);
double Beta,*S=NULL,*Sigma=NULL,kc,*AS=NULL,*BS=NULL;
S=dvector(0,AZNum);Sigma=dvector(0,AZNum);
AS=dvector(0,AZNum);BS=dvector(0,AZNum);
double A=0,B=0;
double **prms1=NULL; prms1=dmatrix(0,HNum,0,obsrvr.FFNum);
double **SPL1=NULL; SPL1=dmatrix(0,HNum,0,obsrvr.FFNum);
enum runtype run_t;
run_t=caseIN.runID;
FILE *fspl1;
double S0,freq;
freq=prop.NB*prop.omega/2/PI;
double **prms2=NULL; prms2=dmatrix(0,HNum,0,obsrvr.FFNum);
double **SPL2=NULL; SPL2=dmatrix(0,HNum,0,obsrvr.FFNum);
FILE *fspl2;
/* */
FILE *fspl,*fsplm,*fsplr;
double **prms=NULL; prms=dmatrix(0,HNum,0,obsrvr.FFNum);
double **SPL=NULL; SPL=dmatrix(0,HNum,0,obsrvr.FFNum);

switch(run_t){
  case(GutinDeming):
    fspl2=fopen("SPLGD.txt","w");
    char fnc[256];strcpy(fnc,caseIN.caseName);
    strcat(fnc,"_SPLGDm.txt");
    fsplm=fopen(fnc,"w");
    for(int m=1;m<HNum+1;m++){
      for(int j=0;j<obsrvr.FFNum;j++){
        Beta = sqrt(1-prop.Mach*prop.Mach);
        S0 = sqrt(pow(obsrvr.FFcoords[j].x[0],2)+Beta*Beta*pow(obsrvr.FFcoords[j].x[1],2));
        kc=m*prop.NB*prop.omega/atm.c0;

        prms2[m][j]=m*prop.NB*prop.omega/(2*PI*sqrt(2)*atm.c0*S0) \
                    *fabs(prop.Thrust*(prop.Mach+obsrvr.FFcoords[j].x[0]/S0)/Beta*Beta \
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
      }
      fprintf(fsplm,"%lf %lf \n",m*freq,SPL2[m][0]);
    }
    fclose(fspl2);
    fclose(fsplm);
    /* second formulation */
    if(DEBUG==1){
      printf("kc, %lf, S0 = %lf\n",kc,S0);
      printf("prms2[0][0] = %lf, prms2[HNum][obsFFNum] = %lf\n",prms2[1][0],prms2[HNum][obsrvr.FFNum-1]);
    }
    break;
  case(GarrickWatkins):
    fspl1=fopen("SPLGW.txt","w");
    char *fn=NULL;strcpy(fn,caseIN.caseName);
    strcat(fn,"_SPLGWm.txt");
    fsplm=fopen(fn,"w");
    for(int m=1;m<HNum+1;m++){
      for(int j=0;j<obsrvr.FFNum;j++){
        /* source position over time */
        for(int k=0;k<AZNum;k++){
          srcPos[k].x[0]=0.0;
          srcPos[k].x[1]=prop.Re*cos(srcTheta[k]);
          srcPos[k].x[2]=prop.Re*sin(srcTheta[k]);
          Beta=sqrt((1-pow(prop.Mach,2)));
          S[k]=sqrt(pow((obsrvr.FFcoords[j].x[0]-srcPos[k].x[0]),2)+
               pow(Beta,2)*(pow((obsrvr.FFcoords[j].x[1]-srcPos[k].x[1]),2)+
               pow((obsrvr.FFcoords[j].x[2]-srcPos[k].x[2]),2)));
            
          Sigma[k]=(prop.Mach*(obsrvr.FFcoords[j].x[0]-srcPos[k].x[0])+S[k])/pow(Beta,2);
          kc=m*prop.NB*prop.omega/atm.c0;
        /* */
        AS[k]=(prop.Thrust*obsrvr.FFcoords[j].x[0]*cos((m+1)*prop.NB*srcTheta[k]+kc*Sigma[k])/pow(S[k],2) \
             +(prop.Thrust*kc*(prop.Mach+obsrvr.FFcoords[j].x[0]/S[k])/pow(Beta,2) \
              -prop.Torque*m*prop.NB/pow(prop.Re,2))*sin(m*prop.NB*srcTheta[k]+kc*Sigma[k]))/S[k];

        BS[k]=(-prop.Thrust*obsrvr.FFcoords[j].x[0]*sin((m+1)*prop.NB*srcTheta[k]+kc*Sigma[k])/pow(S[k],2) \
             +(prop.Thrust*kc*(prop.Mach+obsrvr.FFcoords[j].x[0]/S[k])/pow(Beta,2) \
              -prop.Torque*m*prop.NB/pow(prop.Re,2))*cos(m*prop.NB*srcTheta[k]+kc*Sigma[k]))/S[k];
        }
        A=integTrapz(AS,AZNum); B=integTrapz(BS,AZNum);
        prms1[m][j]=sqrt(2)*sqrt(A*A+B*B)/(8*PI*PI);
        SPL1[m][j]=20*log10(prms1[m][j]/p_ref);
        fprintf(fspl1,"%d %d %lf\n",m,j,SPL1[m][j]);
      }
      fprintf(fsplm,"%lf %lf \n",m*freq,SPL1[m][0]);
    }
    fclose(fspl1);
    fclose(fsplm);

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
    fspl=fopen("SPLBM.txt","w");
    char fnm[256];strcpy(fnm,caseIN.caseName);
    char fnr[256];strcpy(fnr,caseIN.caseName);
    strcat(fnm,"_SPLBMm.txt");
    strcat(fnr,"_SPLBMr.txt");
    fsplm=fopen(fnm,"w"); fsplr=fopen(fnr,"w");
    
    for(int m=1;m<HNum+1;m++){
      for(int j=0;j<obsrvr.FFNum;j++){
        Beta = sqrt(1-prop.Mach*prop.Mach);
        S0 = sqrt(pow(obsrvr.FFcoords[j].x[0],2)+Beta*Beta*pow(obsrvr.FFcoords[j].x[1],2));
        kc=m*prop.NB*prop.omega/atm.c0;
        prop.bladeAngle=prop.bladeAngle*PI/180.0;
        prms[m][j]=1.0/(sqrt(2)*PI*S0)*fabs(prop.Re/(prop.bladeChord*cos(prop.bladeAngle)) \
            *sin(m*prop.NB*prop.bladeChord*cos(prop.bladeAngle)/(2*prop.Re)) \
            *((prop.Mach+obsrvr.FFcoords[j].x[0]/S0)*prop.omega*prop.Thrust/(atm.c0*Beta*Beta) \
            -1/(prop.Re*prop.Re)*prop.Torque)*(besselj(m*prop.NB,kc*prop.Re*obsrvr.FFcoords[j].x[1]/S0) \
            +(pow(Beta,2)*obsrvr.FFcoords[j].x[1]*prop.Re)/(2*S0*S0) \
            *(besselj((m*prop.NB-1),kc*prop.Re*obsrvr.FFcoords[j].x[1]/S0) \
            -besselj((m*prop.NB+1),kc*prop.Re*obsrvr.FFcoords[j].x[1]/S0))));
        SPL[m][j]=20*log10(prms[m][j]/p_ref);
        fprintf(fspl,"%d %d %lf\n",m,j,SPL[m][j]);
        if(m==1) fprintf(fsplr,"%lf %lf\n",(obsrvr.FFcoords[j].x[0]/(PI*5.2)*180.0),SPL[1][j]);
      }
      fprintf(fsplm,"%lf %lf \n",m*freq,SPL[m][0]);
    }
    fclose(fspl);
    fclose(fsplm);
    if(DEBUG==1){
      printf("Sin(30deg) = %lf, sin(pi/6) = %lf, sin(bangle) = %lf\n",sin(30.0),sin(PI/6),sin(prop.bladeAngle*PI/180.0));
    }
    break;
  default:
    APnoise_error("Please specify the run type\n",thisroutine);
}

free_dvector(obsrvr.theta,0,obsrvr.thetaNum);
free_dvector(obsrvr.phi,0,obsrvr.phiNum);
free_dvector(S,0,AZNum);
free_dvector(Sigma,0,AZNum);
free_dvector(AS,0,AZNum);
free_dvector(BS,0,AZNum);
free_d3_tvector(obsrvr.NFcoords,0,obNumAll);
free_d3_tvector(obsrvr.FFcoords,0,FFobsNumAll);
free_d3_tvector(srcPos,0,AZNum);
free_dmatrix(prms1,0,HNum,0,AZNum);
free_dmatrix(SPL1,0,HNum,0,AZNum);
free_dmatrix(prms2,0,HNum,0,AZNum);
free_dmatrix(SPL2,0,HNum,0,AZNum);
free_dmatrix(prms,0,HNum,0,AZNum);
free_dmatrix(SPL,0,HNum,0,AZNum);

/* computing the CPU time */
end=clock();
cpu_time=get_cpu_time(start,end);

return 0;
}


