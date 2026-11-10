
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "structs.h"

#define LIGHTSPEED 3e10

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

