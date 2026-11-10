
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "structs.h"

#define LIGHTSPEED 3e10 
#define PROTON_MASS 1.67e-24
#define COMP_Z 1.0
#define COMP_A 1.0

double get_tgam( struct cell * );
double get_dEbind( double , double , double , double );
double get_diff_dE( struct cell * , struct cell * , double , double );

double get_rho( double , double , double );
double get_temp( double , double );
double get_Press( double , double );
double get_dV( double , double , double );

void timestep( struct domain * theDomain ){

   struct cell * theCells = theDomain->theCells;
   int N = theDomain->N;
   double t = theDomain->t;
   double dt = theDomain->dt;

   int i;
//Nickel Decay
   for( i=0 ; i<N ; ++i ){
      double dm  = theCells[i].dm;
      double dmX = theCells[i].XNi*dm;
      double tgam = get_tgam(theCells+i);

      double dE  = get_dEbind( dmX , t , dt , tgam );

      theCells[i].dE += dE;
   }

//Diffusive Flux
   for( i=0 ; i<N-1 ; ++i ){
      struct cell * cL = theCells+i;
      struct cell * cR = theCells+i+1;
      double dE = get_diff_dE( cL , cR , t , dt );

      theCells[i].dE -= dE;
      theCells[i+1].dE += dE;
   }

//Last Zone
   for(i=N-1;i<N;++i){
      double c = LIGHTSPEED;
      double r = (theCells[i].v+.5*theCells[i].dv)*t;
      double eL = theCells[i  ].eps;

      double F = c*eL;
      double dA = 4.*M_PI*r*r;
      double dE = F*dA*dt;

      theCells[i].dE -= dE;
   }

//Source Term (adiabatic expansion)
   for( i=0 ; i<N ; ++i ){
      double m_p = PROTON_MASS*COMP_A;
      double dMdv = theCells[i].dm/theCells[i].dv;
      double rho = get_rho( dMdv , theCells[i].v , t );
      double n = rho/m_p;

      double eps = theCells[i].eps;
      double T = get_temp( n , eps );
      double P = get_Press( n , T );
      double edot = -3.*P/t;

      double vp = theCells[i].v+.5*theCells[i].dv;
      double vm = theCells[i].v-.5*theCells[i].dv;
      double dV = get_dV( vp, vm , t );
      theCells[i].dE += edot*dt*dV;
   }

//Recompute eps
   for( i=0 ; i<N ; ++i ){
      double vp = theCells[i].v+.5*theCells[i].dv;
      double vm = theCells[i].v-.5*theCells[i].dv;
      double dV = get_dV( vp, vm , t );
      theCells[i].eps = theCells[i].dE/dV;
   }


}
