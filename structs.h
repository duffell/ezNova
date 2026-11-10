
#define NUM_ISOTOPES 8

struct cell{
   double dm;
   double dv;
   double v;
   double eps;
   double dE;
   double tau;
   double XNi;
   double X[NUM_ISOTOPES];
};

struct domain{

   int N;
   int Zs[NUM_ISOTOPES];
   int As[NUM_ISOTOPES];

   double t;
   double t0;
   double dt;

   struct cell * theCells;

};

