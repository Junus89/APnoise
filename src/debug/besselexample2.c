/*****************************
 ******BESSEL FUNCTION********
 ***********SERIES************
 ****************************/
#include<stdio.h>
#include<math.h>
double factorial(int n){
  int i;
  double fact=1;
  for(i=n;i>=1;i--){
    fact=fact*i;
  }
  return fact;
}
int main(){
    FILE *fp=NULL;
    fp=fopen("besselSeriesPlotn0.txt","w");
    double t0,t1,R,sum,x,eps;
    int n;
    printf("Enter the value of n: ");
    scanf("%d",&n);
    printf("Enter the value of x: ");
    scanf("%d",&n);
    printf("Enter the desired accuracy: ");
    scanf("%lf",&eps);
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
        fprintf(fp,"%lf\t%lf\n",x,sum);
  return 0;
}
