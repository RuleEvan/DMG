#include "angular.h"
double triangle(double a, double b, double c) {
  // Computes the triangle coefficients necessary for sixj calculations
  // See Edmonds pg. 99
  double tri = 0.0;
  tri = gsl_sf_gamma(a + b - c + 1.0)*gsl_sf_gamma(a - b + c + 1.0)*gsl_sf_gamma(-a + b + c + 1.0);
  tri *= 1.0/gsl_sf_gamma(a + b + c + 2.0);
  tri = sqrt(tri);

  return tri;
}

double clebsch_gordan(double j1, double j2, double j, double m1, double m2, double m) {
//  printf("Calling CG: %g %g %g %g %g %g\n", j1, m1, j2, m2, j, m);
  // Computes the Clebsch-Gordan coefficients between the uncoupled basis (j1, m1, j2, m2) and
  // the coupled basis (j1,j2; j,m)
  double cg = 0.0;
  if (((m1 + m2) != m) || (j > (j1 + j2)) || (j < abs(j1-j2))) {return cg;}
  double f = 0.0;
  int s_max = MIN(j1 - m1, j - m);
  int s_min = MAX(0, ceil(j -j2 - m1));
//  printf("Min/Max: %d %d\n", s_min, s_max);
  for (int s = s_min; s <= s_max; s++) {
    double d1 = j1 - m1 -s;
    double d2 = j - m - s;
    double d3 = j2 - j + m1 + s;
    if ((d1 < 0) || (d2 < 0) || (d3 < 0)){continue;}
    f += pow(-1.0, s + j1 - m1)*gsl_sf_gamma(j1 + m1 + s +1.0)*gsl_sf_gamma(j2 + j - m1 - s + 1.0)/(gsl_sf_gamma(s + 1.0)*gsl_sf_gamma(d1 + 1.0)*gsl_sf_gamma(d2 + 1.0)*gsl_sf_gamma(d3 + 1.0));
  }
  double w = (2.0*j + 1.0)*gsl_sf_gamma(j1 + j2 - j + 1.0)*gsl_sf_gamma(j1 - m1 + 1.0)*gsl_sf_gamma(j2 - m2 + 1.0)*gsl_sf_gamma(j + m + 1.0)*gsl_sf_gamma(j - m + 1.0)/(gsl_sf_gamma(j1 + j2 + j + 2.0)*gsl_sf_gamma(j1 - j2 + j + 1.0)*gsl_sf_gamma(-j1 + j2 + j + 1.0)*gsl_sf_gamma(j1 + m1 + 1.0)*gsl_sf_gamma(j2 + m2 + 1.0));
  cg = sqrt(w)*f;
//  printf("CG: %g\n", cg);
  return cg;
} 

double three_j(double j1, double j2, double j3, double m1, double m2, double m3) {
  // Computes the Wigner 3J symbol from the corresponding Clebsch-Gordan coefficient
  double three_j = pow(-1.0, j1 - j2 - m3)/sqrt(2.0*j3 + 1.0)*clebsch_gordan(j1, j2, j3, m1, m2, -m3);

  return three_j;
}  

double six_j(double j1, double j2, double j3, double j4, double j5, double j6) {
  // Computes the Wigner six-j symbol using the Racah formual (Edmonds pg. 99)
  double six_j = 0;
  if ((j1 < 0) || (j2 < 0) || (j3 < 0) || (j4 < 0) || (j5 < 0) || (j6 < 0)) {
    printf("Unallowed quantum numbers: %g, %g, %g, %g, %g, %g\n", j1, j2, j3, j4, j5, j6);
    return six_j;
  }
  if ((j1 < abs(j2 - j3)) || (j1 > (j2 + j3)) || (j1 < abs(j5 - j6)) || (j1 > (j5 + j6)) || (j4 < abs(j2 - j6)) || (j4 > (j2 + j6)) || (j4 < abs(j5 - j3)) || (j4 > (j5 + j3))) {return six_j;}
  double t1 = triangle(j1, j2, j3);
  double t2 = triangle(j1, j5, j6);
  double t3 = triangle(j4, j2, j6);
  double t4 = triangle(j4, j5, j3);

  double w = 0.0;
  
  int z;
  for (z = 0; z < 100; z++) {
    double d1 = z - j1 - j2 - j3;
    double d2 = z - j1 - j5 - j6;
    double d3 = z - j4 - j2 - j6;
    double d4 = z - j4 - j5 - j3;
    double d5 = j1 + j2 + j4 + j5 - z;
    double d6 = j2 + j3 + j5 + j6 - z;
    double d7 = j3 + j1 + j6 + j4 - z;
    if ((d1 < 0) || (d2 < 0) || (d3 < 0) || (d4 < 0) || (d5 <0) || (d6 < 0) || (d7 < 0)){continue;}
    double f = gsl_sf_gamma(d1 + 1.0)*gsl_sf_gamma(d2 + 1.0)*gsl_sf_gamma(d3 + 1.0)*gsl_sf_gamma(d4 + 1.0)*gsl_sf_gamma(d5 + 1.0)*gsl_sf_gamma(d6 + 1.0)*gsl_sf_gamma(d7 + 1.0);
    w += pow(-1.0, z)*gsl_sf_gamma(z + 2.0)/f;
  }
  six_j = t1*t2*t3*t4*w;
 
  return six_j;
}

double nine_j(double j11, double j12, double j13, double j21, double j22, double j23, double j31, double j32, double j33) {
  // Computes the Wigner 9J-symbol from the necessary 6J-symbols
  double nine_j = 0.0;
  int i;
 
  for (i = 0; i < 100; i++) {
    double k = i/2.0;
    nine_j += pow(-1.0, i)*(i + 1.0)*six_j(j11, j21, j31, j32, j33, k)*six_j(j12, j22, j32, j21, k, j23)*six_j(j13, j23, j33, k, j11, j12);
  } 
  return nine_j;
}

