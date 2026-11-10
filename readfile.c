
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "structs.h"

#define RAD_CONST  7.5657e-15

void read_header( char * filename , struct domain * theDomain ){

   int N;
   double rin;
   double t;
   int Nc;

   FILE * pFile = fopen(filename,"r");
   fscanf(pFile,"%*s %*s\n%d %lf %lf %d\n", &N , &rin , &t , &Nc );

   printf("#N = %d cells, rin=%.2e cm, t=%.2e days, Nc=%d isotopes\n",N,rin,t/86400.,Nc);
   int Zs[Nc];
   int As[Nc];
   
   int q;
   for( q=0 ; q<Nc ; ++q ){
      fscanf(pFile,"%d.%d ",Zs+q,As+q);
      theDomain->Zs[q] = Zs[q];
      theDomain->As[q] = As[q];
      printf("#Z=%d A=%d\n",Zs[q],As[q]);
   }


   theDomain->N = N;
   theDomain->t0 = t;
   theDomain->t  = t;


   theDomain->theCells = (struct cell *) calloc( N , sizeof(struct cell) );
   struct cell * theCells = theDomain->theCells;

   double rm = rin;
   double vm = rin/t;
   int i;
   for( i=0 ; i<N ; ++i ){
      double rp,vp,rho,T;
      fscanf(pFile, "%lf %lf %lf %lf", &rp , &vp , &rho , &T );
      for( q=0 ; q<Nc ; ++q ){
         fscanf(pFile,"%lf ", theCells[i].X+q );
      }

      double dv = vp-vm;
      double v  = .5*(vp+vm);
      theCells[i].dv = dv;
      theCells[i].v  = v;

      double dm = 4.*M_PI*v*v*dv*t*t*t*rho;
      theCells[i].dm = dm;

      double eps = RAD_CONST*pow(T,4.);
      theCells[i].eps = eps;
      theCells[i].dE  = 4.*M_PI*v*v*dv*t*t*t*eps;

      theCells[i].tau = 0.0;
      theCells[i].XNi = theCells[i].X[Nc-1];

      vm = vp;
      rm = rp;
   }

   fclose(pFile);

}

