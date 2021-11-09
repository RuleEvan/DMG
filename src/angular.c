#include "angular.h"

double clebsch_gordan(double j1, double j2, double j, double m1, double m2, double m) {
/* Computes the Clebsch-Gordan coefficients between the uncoupled basis (j1, m1, j2, m2) and
   the coupled basis (j1,j2; j,m)
*/

  double cg = 0.0;
  if (((m1 + m2) != m) || (j > (j1 + j2)) || (j < abs(j1-j2))) {return cg;}
  cg = pow(-1.0, -j1 + j2 - m)*sqrt(2*j + 1)*three_j(j1, j2, j, m1, m2, -m);  
 
  return cg;
} 

double three_j(double j1, double j2, double j3, double m1, double m2, double m3) {
  // Wrapper to gsl library for 3-j symbol

  double three_j = gsl_sf_coupling_3j((int) 2*j1, (int) 2*j2, (int) 2*j3, (int) 2*m1, (int) 2*m2, (int) 2*m3);
;
  return three_j;
}  

double six_j(double j1, double j2, double j3, double j4, double j5, double j6) {
  // Wrapper to gsl library for 6-j symbol

  double six_j = 0.0;

  six_j = gsl_sf_coupling_6j((int) 2*j1, (int) 2*j2, (int) 2*j3, (int) 2*j4, (int) 2*j5, (int) 2*j6);

  return six_j;
}

double nine_j(double j11, double j12, double j13, double j21, double j22, double j23, double j31, double j32, double j33) {
  // Wrapper to gsl library for 9-j symbol
  double nine_j = 0.0;
   
  nine_j = gsl_sf_coupling_9j((int) 2*j11, (int) 2*j12, (int) 2*j13, (int) 2*j21, (int) 2*j22, (int) 2*j23, (int) 2*j31, (int) 2*j32, (int) 2*j33);

  return nine_j;
}

