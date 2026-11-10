
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "structs.h"


#define LIGHTSPEED 3e10
#define PROTON_MASS 1.67e-24
#define COMP_A 1.0


void read_header( char * , struct domain * );
void report_header( void );
void report( struct domain * );
void timestep( struct domain * );

double get_eps( double , double );
double get_temp( double , double );


double get_rho( double dMdv , double v , double t ){
   return( dMdv/4./M_PI/v/v/t/t/t );
}

double get_dV( double vp , double vm , double t ){
   return( 4./3.*M_PI*(vp*vp*vp-vm*vm*vm)*t*t*t );
}

/*

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


      double T = 2e5*exp(-v/v0)*day/t0;


      theCells[i].dm = dMdv*dv;
      theCells[i].XNi = X;

      double eps = get_eps( 0.0 , T );
      theCells[i].eps = eps;

      double dV = get_dV( vp, vm , t0 );
      theCells[i].dE = theCells[i].eps*dV;
   }

}*/

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

      double dt = dlnt*theDomain.t;
      theDomain.dt = dt;
      ++nstep;

      timestep( &theDomain );

      if(nstep%1000==0) report( &theDomain );
      if(nstep%10000==0) printf("# t=%.2f days\n",theDomain.t/86400.);
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
