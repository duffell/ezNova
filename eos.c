
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define BOLTZMANN 1.38e-16
#define RAD_CONST 7.5657e-15
#define TOL 1e-10


double get_Press( double n , double T ){
   double a = RAD_CONST;
   double k = BOLTZMANN;
   return( n*k*T + a*pow(T,4.)/3. );
}

double get_Prad( double T ){
   double a = RAD_CONST;
   return( a*pow(T,4.)/3. );
}

double get_eps( double n , double T ){
   double k = BOLTZMANN;
   double a = RAD_CONST;

   return( 1.5*n*k*T + a*pow(T,4.) );
}

double dedT( double n , double T ){
   double k = BOLTZMANN;
   double a = RAD_CONST;

   return( 1.5*n*k + 4.*a*pow(T,3.) );
}

double get_temp( double n , double eps ){

   double k = BOLTZMANN;
   double a = RAD_CONST;

   double T1 = eps/(1.5*n*k);
   double T2 = pow( eps/a , .25 );

return(T2);  //Very quick temperature solution if rad-dominated.

   double ns = 2.0;

   double Tg = pow( pow(T1,-ns) + pow(T2,-ns) ,-1./ns );
   double eg = get_eps(n,Tg);
   double de = eg-eps;


   while( fabs(de/eps) > TOL ){
      double deriv = dedT( n , Tg );
      double dT = -de/deriv;
      Tg += dT;
      eg = get_eps(n,Tg);
      de = eg-eps;
   }
   return(Tg);

}

/*
int main(void){
   double Tmin = 1e3;
   double Tmax = 1e9;
   double n = 1e20;

   int N = 100;
   int i;
   for( i=0 ; i<N ; ++i ){
      double x = ((double)i+.5)/(double)N;
      double T = Tmin*pow(Tmax/Tmin,x);
      
      double eps  = get_eps( n , T );
      double Tnew = get_temp( n , eps );
      printf("%e %e %e\n", T , Tnew , eps );
   }

   return(0);
}
*/
