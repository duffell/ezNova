
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "structs.h"

#define TOL 1e-10

#define STEFAN_B   5.67e-5
#define RAD_CONST  7.5657e-15
#define LIGHTSPEED 3e10
#define BOLTZMANN  1.38e-16

#define OPACITY 0.1
#define COMP_Z 1.0
#define COMP_A 1.0
#define PROTON_MASS 1.67e-24

#define BBAND         6.8e14   //Hz
#define BBAND_DNU     1.4e14   //Hz
#define H_OVER_K      4.8e-11  //seconds*Kelvin
#define TP_H_OVER_C2  4.63e-47 //gram*seconds
#define DISTANCE      2e25     //cm

void read_header( char * , struct domain * );
void report_header( void );
void report( struct domain * );
void timestep( struct domain * );

double get_dEbind( double dmX , double t , double dt , double tgam ){

   double day   = 86400;
   double tauNi = 8.8*day;
   double tauCo = 111.3*day;

   double eNi = 9.4e16; //Energy per unit mass of Nickel relative to Iron
   //double eNi = 3.0e16; //Energy per unit mass of Nickel relative to Cobalt
   double eCo = 6.4e16; //Same thing with Cobalt relative to Iron

//Ni -> Co = 1.75 MeV = 3.0e16 erg per unit mass
//Co -> Fe = 3.73 MeV = 6.4e16 erg per unit mass
//Ni -> Fe = 9.4e16 erg per unit mass

   double XNi1 = exp(-t/tauNi);
   double XNi2 = exp(-(t+dt)/tauNi);

   double XCo1 = tauCo/(tauCo-tauNi)*( exp(-t/tauCo) - exp(-t/tauNi) );
   double XCo2 = tauCo/(tauCo-tauNi)*( exp(-(t+dt)/tauCo) - exp(-(t+dt)/tauNi) );

   double dE = dmX*( (XNi1-XNi2)*eNi + (XCo1-XCo2)*eCo );
 
   return( dE*(1.-exp(-tgam*tgam/t/t)) );

}

double get_tgam( struct cell * c ){
   double dm  = c->dm;
   double dv = c->dv;
   double v  = c->v;
   double kappa = 0.03;
   double cs = LIGHTSPEED;

   //double tgam = 0.145*sqrt( dm/dv*kappa*cs/4./M_PI/v/v);
   double tgam = 0.19*sqrt( dm/dv*kappa*cs/4./M_PI/v/v);

   return(tgam);
}

double get_rho( double dMdv , double v , double t ){
   return( dMdv/4./M_PI/v/v/t/t/t );
}

double get_Press( double n , double T ){
   double a = RAD_CONST;
   double k = BOLTZMANN;
   return( 0.0*n*k*T + a*pow(T,4.)/3. );
}

double get_Prad( double T ){
   double a = RAD_CONST;
   return( a*pow(T,4.)/3. );
}

double get_eps( double n , double T ){
   double a = RAD_CONST;
   double k = BOLTZMANN;
   return( 0.0*1.5*n*k*T + a*pow(T,4.) );
}

double dedT( double n , double T ){
   double k = BOLTZMANN;
   double a = RAD_CONST;

   return( 1.5*n*k + 4.*a*pow(T,3.) );
}

double get_temp( double n , double eps ){
   double a = RAD_CONST;
   return( pow(eps/a,.25) );
}

double get_dtau( double dm, double v , double t ){
   double k = OPACITY*COMP_Z/COMP_A;
   return( k*dm/4./M_PI/v/v/t/t );
}

double get_dV( double vp , double vm , double t ){
   return( 4./3.*M_PI*(vp*vp*vp-vm*vm*vm)*t*t*t );
}


double get_diff_dE( struct cell * cL , struct cell * cR , double t , double dt ){
   double c = LIGHTSPEED;
   double m_p = PROTON_MASS*COMP_A;
   double r = (cL->v+.5*cL->dv)*t;

   double dmL = cL->dm;
   double vL  = cL->v;
   double dMdvL = dmL/cL->dv;
   double rhoL = get_rho( dMdvL , vL , t );
   double nL = rhoL/m_p;
   double eL = cL->eps;
   double TL = get_temp( nL , eL );
   double Pr_L = get_Prad( TL );
   double dtauL = .5*get_dtau( dmL , vL , t );

   double dmR = cR->dm;
   double vR  = cR->v;
   double dMdvR = dmR/cR->dv;
   double rhoR = get_rho( dMdvR , vR , t );
   double nR = rhoR/m_p;
   double eR = cR->eps;
   double TR = get_temp( nR , eR );
   double Pr_R = get_Prad( TR );
   double dtauR = .5*get_dtau( dmR , vR , t );

   double dtau = dtauL+dtauR;

   double R = 4.*fabs(log(TR/TL))/dtau;

   double lam = (6. + 3.*R)/(6. + 3.*R + R*R);

   double F = -lam*c*(Pr_R-Pr_L)/dtau;
   double dA = 4.*M_PI*r*r;
   double dE = F*dA*dt;

   return(dE);
}

void initial( struct cell * theCells , int N , double t0 ){

   double c = LIGHTSPEED;
   double Msun = 2e33;
   double day = 86400.;

   double vmax = 3e9;
   double vmin = 1e7;
   double M = 1.3835*Msun;
//   double M = 0.01*Msun;
   double v0 = 2670*1e5; //3.6e8;
//   double v0 = 0.00015*c;

   double T0 = 2.5e7;

   int i;
   for( i=0 ; i<N ; ++i ){
      double xm = ((double)i  )/(double)N;
      double xp = ((double)i+1)/(double)N;
      double vm = vmax*xm;//vmin*pow(vmax/vmin,xm)-vmin;//vmax*xm;
      double vp = vmax*xp;//vmin*pow(vmax/vmin,xp)-vmin;//vmax*xp;
      double dv = vp-vm;
      double v = .5*(vm+vp);
      theCells[i].v  =  v;
      theCells[i].dv = dv;

      double gam = 1./sqrt(1.-v*v/c/c);

      double v1 = 0.0158*c*(0.63)*sqrt(0.075);
      double rho0 = gam*gam*pow(day/t0,3.)*(1./0.075)*sqrt(1./0.63)*6.2e-12/(pow(v/v1,1.0) + pow(v/v1,2.33));

rho0 *= 1.+pow(1e7/v,9.)/(1.+pow(3.2e6/v,8.0));
/* 
//      double dMdv = (M/M_PI)*v*v/v0/v0/v0/( pow(v/v0,2.5) + pow(v/v0,3.5) );
      double dMdv = 4.*M_PI*rho0*v*v*pow(t0,3.);
      double T = T0*pow(v/1e7,0.0)/(1.+pow(v/6e6,7.0))*(day/t0);
//      double T = T0*pow(v/c,-7./8.);
      double X = 0.0;
      //if( v<6e6 ) X = 1000.0;
*/

      double dMdv = .5*M*v*v/v0/v0/v0*exp(-v/v0);
      double vrho = 7e7;
      if( v<vrho ) dMdv = .5*M*v*v/v0/v0/v0*exp(-vrho/v0);

      double delta = 0.0;
      if(v>2e8) delta = -0.075;
      if(v>9.9e8) delta += 0.3;
      dMdv*= 1. + delta;


      double X = 0.0;
      //if( v<2.95*v0 ) X = 1.0;
      if( v<9.9e8 ) X = 1.0;
      double vin = 3.5e8;
      if( v<vin ) X *= exp((v-vin)/1.7e8);
      double Xin = v/4.5e8;
      if( X > Xin ) X = Xin;

//      double m_p = PROTON_MASS*COMP_A;

      double T = 2e5*exp(-v/v0)*day/t0;


      theCells[i].dm = dMdv*dv;
      theCells[i].XNi = X;

      double eps = get_eps( 0.0 , T );
      theCells[i].eps = eps;

      double dV = get_dV( vp, vm , t0 );
      theCells[i].dE = theCells[i].eps*dV;
   }

}

int main( int argc , char * argv[] ){

   if( argc<2 ){ printf("Please specify the input file name.\n");  return(0);}

   report_header();

   struct domain theDomain;
   read_header(argv[1], &theDomain );

   int N = 100;
   N = theDomain.N;
   double vmax = 3e9;
   vmax = theDomain.theCells[N-1].v;

   double CFL = 0.003;
   double dlnt = CFL*vmax/(double)N/LIGHTSPEED;


   double day = 86400;
   double t0 = 2.0*day;
   t0 = theDomain.t0;
   double tmax = 150.*day;

 
   struct cell * theCells = theDomain.theCells;


//   initial( theCells , N , t0 );


   double MNi = 0.0;
   double MTot = 0.0;
   double ETot = 0.0;
   double Eth  = 0.0;
   int i;
   for( i=0; i<N ; ++i ){
      double v = theCells[i].v;
      Eth  += theCells[i].dE;
      MTot += theCells[i].dm;
      ETot += .5*theCells[i].dm*v*v;
      MNi += theCells[i].dm*theCells[i].XNi;
   }
   printf("#Kinetic Energy = %.4e erg\n",ETot);
   printf("#Thermal Energy at 1 day = %.4e erg\n",Eth*(t0/day));
   printf("#Total Mass = %.4e Msun\n",MTot/2e33);
   printf("#Nickel Mass = %.4e Msun\n",MNi/2e33);

   int nstep = 0;   


   while( theDomain.t<tmax ){

      double t = theDomain.t;
      double dt = dlnt*t;
      theDomain.dt = dt;
      ++nstep;

      timestep( &theDomain );

//Report
      if(nstep%1000==0) report( &theDomain );
      theDomain.t += dt;

   }



   FILE *pFile = fopen("output.dat","w");

   for( i=0 ; i<N ; ++i ){
      double v = theDomain.theCells[i].v;
      double r = v*theDomain.t;
      double dMdv = theDomain.theCells[i].dm/theDomain.theCells[i].dv;
      double rho = get_rho(dMdv,v, theDomain.t);
      double m_p = PROTON_MASS*COMP_A;
      double n = rho/m_p;
      double T = get_temp( n , theCells[i].eps );
      double tau = theDomain.theCells[i].tau;
      double X = theDomain.theCells[i].XNi;
      fprintf(pFile,"%e %e %e %e %e\n",v,rho,T,tau,X); 
   }


   free(theDomain.theCells);
   fclose(pFile);
}
