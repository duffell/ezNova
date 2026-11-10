
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

double me   = 9.1093837e-28;
double c    = 2.99792e10;
double hbar = 1.05457182e-27;
double k    = 1.380649e-16;
double a    = 7.5657e-15;

double get_P_NR( double n ){
   return( 0.2*pow(3.*M_PI*M_PI,2./3.)*hbar*hbar/me*pow(n,5./3.) );
}

double get_P_UR( double n ){
   return( 0.25*pow(3.*M_PI*M_PI,1./3.)*hbar*c*pow(n,4./3.) );
}

double get_P_supercrude( double n ){
   double n_0 = 125./64.*pow(me*c,3.)/3./M_PI/M_PI/hbar/hbar/hbar;
   double P_0 = 625./256.*c/12.*pow(me*c,4.)/M_PI/M_PI/hbar/hbar/hbar;

   return( P_0/sqrt( pow(n/n_0,-10./3.) + pow(n/n_0,-8./3.) ) );
}

double get_eps_NR( double n ){
   return( 0.3*pow(3.*M_PI*M_PI,2./3.)*hbar*hbar/me*pow(n,5./3.) );
}

double get_eps_UR( double n ){
   return( 0.75*pow(3.*M_PI*M_PI,1./3.)*hbar*c*pow(n,4./3.) );
}

double get_eps_supercrude( double n ){
   double n_1 = 125./8.*1./3./M_PI/M_PI*pow(me*c/hbar,3.);
   double e_1 = 625./64./M_PI/M_PI*pow(me*c/hbar,3.)*me*c*c;
   
   return( e_1/sqrt( pow( n/n_1 , -10./3. ) + pow(n/n_1, -8./3.) ) );
}

double get_P_chandra( double n ){
   double n_R = 5.8651566e29;
   double P_R = 1.44055972e24;
   double x = pow(n/n_R,1./3.);

   double f = x*x*x*sqrt(1.+x*x)/3. - x*(2.*x*x+1.)*sqrt(1.+x*x)/8. + asinh(x)/8.;

   return(P_R*f);
}

double get_eps_chandra( double n ){
   double n_R = 5.8651566e29;
   double P_R = 1.44055972e24;
   double x = pow(n/n_R,1./3.);
   
   double f = x*(2.*x*x+1.)*sqrt(1.+x*x)/8. - asinh(x)/8. - x*x*x/3.;

   return(P_R*f);
}

double get_P_Therm( double n , double T ){
   return( n*k*T + a*T*T*T*T/3. );
}

double get_eps_Therm( double n , double T){
   return( 1.5*n*k*T + a*T*T*T*T );
}

double superCrude_get_P( double n , double T ){
   double P_Th = get_P_Therm( n , T );
   double P_d  = get_P_chandra( n );

//   return( sqrt( P_Th*P_Th + P_d*P_d ) );
   return( P_Th + P_d );
}

double superCrude_get_eps( double n , double T ){
   double eps_Th = get_eps_Therm( n , T );
   double eps_d  = get_eps_chandra( n );
   
   return( eps_Th + eps_d );
}

double get_Phi( double n , double T ){
   double phi_d  = 0.6*pow(3.*M_PI*M_PI*n,2./3.)*hbar*hbar/me;
   double phi_th = -((1./3.)*a/n*T*T*T*T + 1.5*k*T*log(T)-k*T*log(n));

   return( phi_d + phi_th );
}

int main(void){

   double n_min = 1e20;
   double n_max = 1e35;

   double T = 1e8;

   int N = 1000;
   double dx = 1./(double)N;

   int i;
   for( i=0 ; i<N ; ++i ){
      double x = ((double)i+0.5)*dx;
      double n = n_min*pow(n_max/n_min,x);
      
//      double P_NR = get_eps_NR( n );
//      double P_UR = get_eps_UR( n );
//      double P_SC = get_eps_supercrude( n );
//      double P_CH = get_eps_chandra( n );
      double phi = get_Phi( n , T );
      double dphi = get_Phi( n*(1.0001) , T ) - phi;
      double P_phi = n*dphi/0.0001;
      
      double P_SC = superCrude_get_P( n , T );
      double eps_SC = superCrude_get_eps( n , T );

      printf("%e %e %e\n",n,P_SC,eps_SC);
   }

   return(0);
}
