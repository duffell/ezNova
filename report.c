
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "structs.h"

#define STEFAN_B   5.67e-5
#define RAD_CONST  7.5657e-15
#define LIGHTSPEED 3e10
#define BOLTZMANN  1.38e-16

#define BBAND         6.8e14   //Hz
#define BBAND_DNU     1.4e14   //Hz
#define H_OVER_K      4.8e-11  //seconds*Kelvin
#define TP_H_OVER_C2  4.63e-47 //gram*seconds
#define DISTANCE      2e25     //cm


double get_dtau( double , double , double );
//void find_phot( struct cell * , int , double , int * , double * );

double get_diff_dE( struct cell * , struct cell * , double , double );
double get_tgam( struct cell * );
double get_dEbind( double , double , double , double );
double get_temp( double , double );


void fill_taus( struct cell * theCells, int N , double t ){
   int i;
   double tau = 0.0;
   for( i=N-1; i>=0 ; --i ){
      struct cell * c = theCells+i;
      double dm = c->dm;
      double v  = c->v;
      double dtau = get_dtau( dm , v , t );
      tau += .5*dtau;
      c->tau = tau;
      tau += .5*dtau;
   }
}

void find_phot( struct cell * theCells , int N , double t , int * i , double * Tp ){
   
   int j=N-2;
   double c = LIGHTSPEED;
   double tau = theCells[j].tau;
   double tphot = 1./3.;
   while(tau<tphot && j>0){
      --j;
      tau = theCells[j].tau;
   }
   *i = j;
   double eps = theCells[j].eps;
   double v = theCells[j].v;
  
   *Tp = get_temp( 1.0 , eps );//pow(c/v,0.25);

}


void report_header( void ){

   FILE * pFile = fopen("report.dat","w");
   fprintf(pFile,"#1:time    2:i_phot    3:Rphot    4:Tphot    5:L    6:LNi    7:Lnu    8:LNi_nu\n");
   fclose(pFile);

}

void report( struct domain * theDomain ){

   double day = 86400.;
   struct cell * theCells = theDomain->theCells;
   int N = theDomain->N;
   double t = theDomain->t;
   double dt = theDomain->dt;

   fill_taus( theCells , N , t );
   int iphot;
   double Tphot;
   find_phot( theCells , N , t , &iphot , &Tphot );


   double sig = STEFAN_B;
      //double L = 4.*M_PI*Rphot*Rphot*sig*pow(T,4.);

   struct cell * cL = theCells+iphot;
   struct cell * cR = theCells+iphot+1;
   double L = get_diff_dE( cL , cR , t , dt )/dt;

   double Rphot = (cL->v + .5*cL->dv)*t;


   double nu = BBAND;
   double LNi = 0.0;
   double LNi_nu = 0.0;

   int i;
   for( i=iphot+1 ; i<N ; ++i ){
      double dmX = theCells[i].dm*theCells[i].XNi;
      double tgam = get_tgam(theCells+i);
      double dE = get_dEbind( dmX , t , dt , tgam );
      LNi += dE/dt;

      double eps = theCells[i].eps;
      double T = get_temp( 1.0 , eps );
      double u = H_OVER_K*nu/T;
      double Lnu_over_L = 15./pow(M_PI,4.)*H_OVER_K/T*u*u*u/(exp(u)-1.);

      LNi_nu += dE/dt*Lnu_over_L*BBAND_DNU;
   }


      //double C = TP_H_OVER_C2;
   double u = H_OVER_K*nu/Tphot;
   double d = DISTANCE;


   double Lnu_over_L = 15./pow(M_PI,4.)*H_OVER_K/Tphot*u*u*u/(exp(u)-1.);
   double Lnu = L*Lnu_over_L*BBAND_DNU;
   L += LNi;

   FILE * pFile = fopen("report.dat","a");
   fprintf(pFile,"%e %d %e %e %e %e %e %e\n",t/day,iphot,Rphot,Tphot,L,LNi,Lnu,LNi_nu);
   fclose(pFile);
}
