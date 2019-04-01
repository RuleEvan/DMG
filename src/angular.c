#include "angular.h"

double clebsch_gordan(double j1, double j2, double j, double m1, double m2, double m) {
/* Computes the Clebsch-Gordan coefficients between the uncoupled basis (j1, m1, j2, m2) and
   the coupled basis (j1,j2; j,m)
*/

  double cg = 0.0;
  if (((m1 + m2) != m) || (j > (j1 + j2)) || (j < abs(j1-j2))) {return cg;}
  cg = pow(-1.0, -j1 + j2 - m)*sqrt(2*j + 1)*three_j(j1, j2, j, m1, m2, -m);  
//  printf("%g, %g, %g, %g, %g, %g, %g\n", j1, m1, j2, m2, j, m, cg);
  return cg;
} 

double three_j(double j1, double j2, double j3, double m1, double m2, double m3) {
/* Computes the Wigner 3J symbol from the corresponding Clebsch-Gordan coefficient
*/

  double three_j = gsl_sf_coupling_3j((int) 2*j1, (int) 2*j2, (int) 2*j3, (int) 2*m1, (int) 2*m2, (int) 2*m3);
;
  return three_j;
}  

double six_j(double j1, double j2, double j3, double j4, double j5, double j6) {
  double six_j = 0.0;

  six_j = gsl_sf_coupling_6j((int) 2*j1, (int) 2*j2, (int) 2*j3, (int) 2*j4, (int) 2*j5, (int) 2*j6);

  return six_j;
}

double nine_j(double j11, double j12, double j13, double j21, double j22, double j23, double j31, double j32, double j33) {
  // Computes the Wigner 9J-symbol from the necessary 6J-symbols
  double nine_j = 0.0;
   
  nine_j = gsl_sf_coupling_9j((int) 2*j11, (int) 2*j12, (int) 2*j13, (int) 2*j21, (int) 2*j22, (int) 2*j23, (int) 2*j31, (int) 2*j32, (int) 2*j33);
//  printf("%g %g %g %g %g %g %g %g %g %g\n", j11, j12, j13, j21, j22, j23, j31, j32, j33, nine_j); 
  return nine_j;
}

/* Deprecated code

double six_j(double j1, double j2, double j3, double j4, double j5, double j6) {
 Computes the Wigner six-j symbol using the Racah formual (Edmonds pg. 99)


  double six_j = 0;
  if ((j1 < 0) || (j2 < 0) || (j3 < 0) || (j4 < 0) || (j5 < 0) || (j6 < 0)) {
    printf("Unallowed quantum numbers: %g, %g, %g, %g, %g, %g\n", j1, j2, j3, j4, j5, j6);
    return six_j;
  }
  if ((j1 < fabs(j2 - j3)) || (j1 > (j2 + j3)) || (j1 < fabs(j5 - j6)) || (j1 > (j5 + j6)) || (j4 < fabs(j2 - j6)) || (j4 > (j2 + j6)) || (j4 < fabs(j5 - j3)) || (j4 > (j5 + j3))) {return six_j;}
  double t1 = triangle(j1, j2, j3);
  double t2 = triangle(j1, j5, j6);
  double t3 = triangle(j4, j2, j6);
  double t4 = triangle(j4, j5, j3);

  double w = 0.0;
  
  int z_min = MAX(1, j1 + j2 + j3);
  z_min = MAX(z_min, j1 + j5 + j6);
  z_min = MAX(z_min, j4 + j2 + j6);
  z_min = MAX(z_min, j4 + j5 + j3);
 
  int z_max = MIN(j1 + j2 + j4 + j5, j2 + j3 + j5 + j6);
  z_max = MIN(z_max, j3 + j1 + j6 + j4);
  if (z_max < z_min) {return 0.0;}

  for (int z = z_min; z <= z_max; z++) {
    double d1 = z - j1 - j2 - j3;
    double d2 = z - j1 - j5 - j6;
    double d3 = z - j4 - j2 - j6;
    double d4 = z - j4 - j5 - j3;
    double d5 = j1 + j2 + j4 + j5 - z;
    double d6 = j2 + j3 + j5 + j6 - z;
    double d7 = j3 + j1 + j6 + j4 - z;
    double f = gsl_sf_gamma(d1 + 1.0)*gsl_sf_gamma(d2 + 1.0)*gsl_sf_gamma(d3 + 1.0)*gsl_sf_gamma(d4 + 1.0)*gsl_sf_gamma(d5 + 1.0)*gsl_sf_gamma(d6 + 1.0)*gsl_sf_gamma(d7 + 1.0);
    w += pow(-1.0, z)*gsl_sf_gamma(z + 2.0)/f;
  }
  six_j = t1*t2*t3*t4*w;
//  printf("Six_j: %g, %g, %g, %g, %g, %g, %g\n", j1, j2, j3, j4, j5, j6, six_j);
  return six_j;
}

double clebsch_gordan(double j1, double j2, double j, double m1, double m2, double m) {
 Computes the Clebsch-Gordan coefficients between the uncoupled basis (j1, m1, j2, m2) and
  the coupled basis (j1,j2; j,m)


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
  return cg;
} 

double triangle(double a, double b, double c) {
 Computes the triangle coefficients necessary for sixj calculations
   See Edmonds pg. 99


  double tri = 0.0;
  tri = gsl_sf_gamma(a + b - c + 1.0)*gsl_sf_gamma(a - b + c + 1.0)*gsl_sf_gamma(-a + b + c + 1.0);
  tri *= 1.0/gsl_sf_gamma(a + b + c + 2.0);
  tri = sqrt(tri);

  return tri;
}

*/
