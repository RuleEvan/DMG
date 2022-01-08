#include "angular.h"

double log_factorial(unsigned int n) {
  if (n < 0) {printf("Error: Log Factorial expects non-negative argument\n"); exit(0);}
  if (n == 0 || n == 1) {return 0;}
 
  double logf = 0;

  for (int i = 2; i <= n; i++) {
    logf += log(i);
  }

  return logf;
}

int test_log_factorial() {
  int pass = 1;
  if (log_factorial(0) != 0) {pass = 0;}
  if (log_factorial(1) != 0) {pass = 0;}
  if (log_factorial(2) != log(2)) {pass = 0;}
  if (log_factorial(17) != log(355687428096000)) {pass = 0;}
  
  if (pass == 1) {
    printf("Test Log Factorial: Pass\n");
  } else {
    printf("Test Log Factorial: Fail\n");
  }

  return pass;
}

int test_clebsch_gordan() {
  int pass = 1;
  if (clebsch_gordan(0.0,0.0,0.0,0.0,0.0,0.0) != 1) {pass = 0;}
  if (fabs(clebsch_gordan(0.5,0.5,1.0,0.5,-0.5,0.0) - 1.0/sqrt(2.0) > pow(10, -6))) {pass = 0;}
  if (fabs(clebsch_gordan(0.5,2.0,1.5,0.5,-1.0,-0.5) - sqrt(3.0/5.0)) > pow(10, -6)) {pass = 0;} 
  if (clebsch_gordan(0.5,2.0,1.5,0.5,1.0,-0.5) != 0.0) {pass = 0;}
  if (clebsch_gordan(2.0,2.0,6.0,1.0,-1.0,0.0) != 0.0) {pass = 0;}
  if (fabs(clebsch_gordan(13.5,12.5,24.0,12.5,-11.5,1.0) - 191/(10*sqrt(281132955186))) > pow(10, -6)) {pass = 0;}
  
  if (pass == 1) {
    printf("Test Clebsch-Gordan: Pass\n");
  } else {
    printf("Test Clebsch-Gordan: Fail\n");
  }

  return pass;
}


double clebsch_gordan(double j1, double j2, double j, double m1, double m2, double m) {
/* Computes the Clebsch-Gordan coefficients between the uncoupled basis (j1, m1, j2, m2) and
   the coupled basis (j1,j2; j,m)
*/

  double cg = 0.0;
  if (j1 < 0 || j2 < 0 || j < 0) {printf("Error negative angular momentum in Clebsch-Gordan: j1 = %g j2 = %g j = %g\n", j1, j2, j); exit(0);}
  if ((int) (2*j1) % 2 != (int) (2*fabs(m1)) % 2) {printf("Error: unphysical angular momentum j1 = %g m1 = %g\n", j1, m1); exit(0);}
  if ((int) (2*j2) % 2 != (int) (2*fabs(m2)) % 2) {printf("Error: unphysical angular momentum j2 = %g m2 = %g\n", j2, m2); exit(0);}
  if ((int) (2*j) % 2 != (int) (2*fabs(m)) % 2) {printf("Error: unphysical angular momentum j = %g m = %g\n", j, m); exit(0);}

  if (((m1 + m2) != m) || (j > (j1 + j2)) || (j < abs(j1-j2))) {return cg;}
  
  int u1 = (int) (j1 + j2 - j);
  int u2 = (int) (j1 - j2 + j);
  int u3 = (int) (-j1 + j2 + j);
  int u4 = (int) (j1 + j + j2 + 1);


  int w1 = (int) (j1 + m1);
  int w2 = (int) (j1 - m1);
  int w3 = (int) (j2 + m2);
  int w4 = (int) (j2 - m2);
  int w5 = (int) (j + m);
  int w6 = (int) (j - m);


  double delta = sqrt((2*j + 1))*exp(0.5*(log_factorial(u1) + log_factorial(w2) + log_factorial(w4) + log_factorial(w5) + log_factorial(w6) - log_factorial(u4) - log_factorial(u2) - log_factorial(u3) - log_factorial(w1) - log_factorial(w3)));
  int k1 = (int) (j2 + j - m1);
  int k2 = (int) (j2 - j + m1);

  int smin = MAX(0, -k2);
  int smax = MIN(w2, w6);

  for (int s = smin; s <= smax; s++) {
    cg += pow(-1.0, s + w2)*exp(log_factorial(w1 + s) + log_factorial(k1 - s) - log_factorial(s) - log_factorial(w2 - s) - log_factorial(w6 - s) - log_factorial(k2 + s));

  }

  cg *= delta;

  return cg;
} 

