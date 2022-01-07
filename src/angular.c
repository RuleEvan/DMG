#include "angular.h"

double factorial(int n) {
  if (n == 0 || n == 1) {return 1;}

  return n*factorial(n - 1);
}

double clebsch_gordan(double j1, double j2, double j, double m1, double m2, double m) {
/* Computes the Clebsch-Gordan coefficients between the uncoupled basis (j1, m1, j2, m2) and
   the coupled basis (j1,j2; j,m)
*/

  double cg = 0.0;
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


  double delta = sqrt((2*j + 1)*factorial(u1)*factorial(w2)*factorial(w4)*factorial(w5)*factorial(w6)/(factorial(u4)*factorial(u2)*factorial(u3)*factorial(w1)*factorial(w3)));
  int k1 = (int) (j2 + j - m1);
  int k2 = (int) (j2 - j + m1);

  int smin = MAX(0, -k2);
  int smax = MIN(w2, w6);

  for (int s = smin; s <= smax; s++) {
    cg += pow(-1.0, s + w2)*factorial(w1 + s)*factorial(k1 - s)/(factorial(s)*factorial(w2 - s)*factorial(w6 - s)*factorial(k2 + s));

  }

  cg *= delta;

  return cg;
} 

