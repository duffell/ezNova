
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "structs.h"

#define LIGHTSPEED 3e10
#define PROTON_MASS 1.67e-24
   
#define COMP_A 1.0

double get_rho( double , double , double );
double get_temp( double , double );
double get_Prad( double );

double get_dtau( double dm, double v , double t , double kappa ){
   double k = kappa;
   return( k*dm/4./M_PI/v/v/t/t );
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
   double dtauL = .5*get_dtau( dmL , vL , t , cL->kappa );

   double dmR = cR->dm;
   double vR  = cR->v;
   double dMdvR = dmR/cR->dv;
   double rhoR = get_rho( dMdvR , vR , t );
   double nR = rhoR/m_p;
   double eR = cR->eps;
   double TR = get_temp( nR , eR );
   double Pr_R = get_Prad( TR );
   double dtauR = .5*get_dtau( dmR , vR , t , cR->kappa );

   double dtau = dtauL+dtauR;

   double R = 4.*fabs(log(TR/TL))/dtau;

   double lam = (6. + 3.*R)/(6. + 3.*R + R*R);

   double F = -lam*c*(Pr_R-Pr_L)/dtau;
   double dA = 4.*M_PI*r*r;
   double dE = F*dA*dt;

   return(dE);
}
