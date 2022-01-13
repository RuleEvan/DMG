#include "slater.h"

unsigned int p_step(int n_s, int n_p, int *m_p) {
/* Compute the p-coefficient (see Whitehead) associated with a
   given Slater determinant

  Inputs(s):
    int n_s: number of single-particle states
    int n_p: number of particles
    int array m_p: for each particle, specifies the m-value 
                   of the orbital that it occupies
  Output(s): 
    int p: p-coefficient of the given Slater determinant
*/
  unsigned int p = n_choose_k(n_s, n_p);
  for (int i = 0; i < n_p; i++) {
    p -= n_choose_k(n_s-m_p[i], n_p - i);
  }
  return p;
}

int test_slater_routines() {
  int pass = 1;

  int np = 5;
  int ns = 10;
  if (get_num_sds(ns, np) != 252) {pass = 0;}

  int *occ = malloc(np*sizeof(int));
  int *lshell = malloc(ns*sizeof(int));
  int *wshell = malloc(ns*sizeof(int));

  int *jzshell = malloc(ns*sizeof(int));

  lshell[0] = 1;
  lshell[1] = 1;
  lshell[2] = 1;
  lshell[3] = 1;
  lshell[4] = 2;
  lshell[5] = 2;
  lshell[6] = 2;
  lshell[7] = 2;
  lshell[8] = 2;
  lshell[9] = 2;

  jzshell[0] = -3;
  jzshell[1] = -1;
  jzshell[2] = 1;
  jzshell[3] = 3;
  jzshell[4] = -5;
  jzshell[5] = -3;
  jzshell[6] = -1;
  jzshell[7] = 1;
  jzshell[8] = 3;
  jzshell[9] = 5;

  wshell[0] = 1;
  wshell[1] = 1;
  wshell[2] = 1;
  wshell[3] = 1;
  wshell[4] = 3;
  wshell[5] = 3;
  wshell[6] = 3;
  wshell[7] = 3;
  wshell[8] = 3;
  wshell[9] = 3;

  float minmj = min_mj(ns, np, jzshell);
  if (minmj != -6.5) {pass = 0;}
  float maxmj = max_mj(ns, np, jzshell);
  if (maxmj != 6.5) {pass = 0;}
  
  int wmin = min_w(ns, np, wshell);
  if (wmin != 7) {pass = 0;} 


  occ[0] = 1;
  occ[1] = 2;
  occ[2] = 3;
  occ[3] = 4;
  occ[4] = 5;

  unsigned int pi = p_step(ns, np, occ);
  unsigned int pf;
  int phase = 1;
  int j_min =  j_min_from_p(ns, np, pi);
  if (j_min != 1) {pass = 0;}

  int par = parity_from_p(pi, ns, np, lshell);
  if (par != 1) {pass = 0;}

  int w = w_from_p(pi, ns, np, wshell);
  if (w != 7) {pass = 0;}

  float mj = m_from_p(pi, ns, np, jzshell);
  if (mj != -2.5) {pass = 0;}


  if (pi != 1) {pass = 0;};
  for (int n = 1; n <= np; n++) {
    pf = a_op_dag(ns, np, pi, n, &phase, j_min);
    if (pf != 0) {pass = 0;}
  }
  phase = 1;
  pf = a_op_dag(ns, np, pi, np + 1, &phase, j_min);
  if (pf != 1) {pass = 0;}
  if (phase != -1) {pass = 0;}
  phase = 1;
  pf = a_op_dag(ns, np, pi, np + 2, &phase, j_min);
  if (pf != 2) {pass = 0;}
  if (phase != -1) {pass = 0;}

  for (int n = np + 1; n <= ns; n++) {
    pf = a_op(ns, np, pi, n, &phase, j_min);
    if (pf != 0) {pass = 0;}
  }
  phase = 1;
  pf = a_op(ns, np, pi, np, &phase, j_min);
  if (pf != 1) {pass = 0;}
  if (phase != 1) {pass = 0;}
  par = parity_from_p(pf, ns, np - 1, lshell);
  if (par != 1) {pass = 0;}
  w = w_from_p(pf, ns, np - 1, wshell);
  if (w != 4) {pass = 0;}
  mj = m_from_p(pf, ns, np - 1, jzshell);
  if (mj != 0) {pass = 0;}


  phase = 1;
  pf = a_op(ns, np, pi, np - 1, &phase, j_min);
  if (pf != 2) {pass = 0;}
  if (phase != -1) {pass = 0;}
  par = parity_from_p(pf, ns, np - 1, lshell);
  if (par != -1) {pass = 0;}
  w = w_from_p(pf, ns, np - 1, wshell);
  if (w != 6) {pass = 0;}
  mj = m_from_p(pf, ns, np - 1, jzshell);
  if (mj != -4) {pass = 0;}


  occ[0] = 6;
  occ[1] = 7;
  occ[2] = 8;
  occ[3] = 9;
  occ[4] = 10;


  pi = p_step(ns, np, occ);
  if (pi != 252) {pass = 0;};
  par = parity_from_p(pi, ns, np, lshell);
  if (par != 1) {pass = 0;}
  w = w_from_p(pi, ns, np, wshell);
  if (w != 15) {pass = 0;}
  mj = m_from_p(pi, ns, np, jzshell);
  if (mj != 2.5) {pass = 0;}

  occ[0] = 2;
  occ[1] = 4;
  occ[2] = 5;
  occ[3] = 7;
  occ[4] = 10;

  pi = p_step(ns, np, occ);
  if (pi != 168) {pass = 0;}
  par = parity_from_p(pi, ns, np, lshell);
  if (par != 1) {pass = 0;}
  w = w_from_p(pi, ns, np, wshell);
  if (w != 11) {pass = 0;}
  mj = m_from_p(pi, ns, np, jzshell);
  if (mj != 0.5) {pass = 0;}

  if (pass == 1) {
    printf("Test Slater Routines: Pass\n");
  } else {
    printf("Test Slater Routines: Fail\n");
  }

  return pass;
}



unsigned int a_op_dag(int n_s, int n_p, unsigned int p, int n_op, int *phase, int j_min) {
/* Acts an annihilation operator on the orbital in position n_op
   Uses a similar algorithm to orbitals_from_p to reconstruct the occupied orbitals from the p-coefficient

  Input(s):
    int n_s: number of single-particle states
    int n_p: number of particles
    unsigned int p: initial SD p-coefficient
    int n_op: orbital on which the operator acts, labeling starts at 1 in accordance with Whitehead
    int *phase: phase +/- 1 that results from anticommutations
    int j_min: first occupied orbital, can be used to save time when it is known that first j_min orbitals are empty

  Output(s):
    unsigned int pf: resulting SD p-coefficient
*/
  *phase = 1;
  // We use the method of Whitehead to mod out by binomial coefficients to find occupied orbitals
  unsigned int q = n_choose_k(n_s, n_p) - p;

  // k tracks the number of particles found in the initial SD
  int ki = 0;
  // kf tracks the number of particles in the final SD
  int kf = 0;

  // Initialize the final state p-coefficient
  unsigned int pf = n_choose_k(n_s, n_p + 1);
  for (int j = 0; j <= n_s; j++) {
    // Check if all of the particles have been found
    // This inequality is satisfied if the j-th orbital is occupied
    if ((q >= n_choose_k(n_s - j, n_p - ki)) && (q < n_choose_k(n_s - (j - 1), n_p - ki))) {
      if (j != n_op) { // Found orbital to keep in the final SD
        pf -= n_choose_k(n_s - j, n_p - kf + 1);
        kf++;
      } else { // State is already occupied
        return 0;
      }
      if (j < n_op) { // We must anticommute the annihilation operator past this creation operator
        *phase *= -1;
      }
      q -= n_choose_k(n_s - j, n_p  - ki); 
      ki++;
    } else if (j == n_op) { // Target state is empty so we add it in
      pf -= n_choose_k(n_s - j, n_p - kf + 1);
      kf++;
    }
  }
  return pf;
}


unsigned int a_op(int n_s, int n_p, unsigned int p, int n_op, int *phase, int j_min) {
  if (n_p == 0) {return 0;}
/* Acts an annihilation operator on the orbital in position n_op
   Uses a similar algorithm to orbitals_from_p to reconstruct the occupied orbitals from the p-coefficient

  Input(s):
    int n_s: number of single-particle states
    int n_p: number of particles
    unsigned int p: initial SD p-coefficient
    int n_op: orbital on which the operator acts, labeling starts at 1 in accordance with Whitehead
    int *phase: phase +/- 1 that results from anticommutations
    int j_min: first occupied orbital, can be used to save time when it is known that first j_min orbitals are empty

  Output(s):
    unsigned int pf: resulting SD p-coefficient
*/
  *phase = 1;
  // We use the method of Whitehead to mod out by binomial coefficients to find occupied orbitals
  unsigned int q = n_choose_k(n_s, n_p) - p;

  // k tracks the number of particles found in the initial SD
  int ki = 0;
  // kf tracks the number of particles in the final SD
  int kf = 0;

  // Initialize the final state p-coefficient
  unsigned int pf = n_choose_k(n_s, n_p - 1);
  for (int j = j_min; j <= n_s; j++) {
    // Check if all of the particles have been found
    if (ki == n_p) {
      // If all initial state particles are found before the position of the annihilation operator at n_op,
      // then we know that this state is unoccupied and we can return zero
      if (j <= n_op) {return 0;} else {break;}
    }
    // This inequality is satisfied if the j-th orbital is occupied
    if ((q >= n_choose_k(n_s - j, n_p - ki)) && (q < n_choose_k(n_s - (j - 1), n_p - ki))) {
      if (j != n_op) { // Found orbital to keep in the final SD
        pf -= n_choose_k(n_s - j, n_p - kf - 1);
        kf++;
      }
      if (j < n_op) { // We must anticommute the annihilation operator past this creation operator
        *phase *= -1;
      }
      q -= n_choose_k(n_s - j, n_p  - ki); 
      ki++;
    } else if (j == n_op) { // Target state is empty so operator annihilates the initial SD
      return 0;
    }
  }
  return pf;
}

unsigned int get_num_sds(int n_s, int n_p) {
  /* Given the number of single-particle states and number of particles, computes the number of possible SDs
     
  Input(s):
    int n_s: number of single-particle states
    int n_p: number of particles

  Output(s):
    unsigned int n_sds: number of Slater determinants
*/

  unsigned int n_sds = 1;
  if (n_p <= 0) {return 1;}
  int* max_state = (int*) malloc(sizeof(int)*n_p);
  for (int i = 0; i < n_p; i++) {
    max_state[i] = n_s - (n_p - i - 1);
  }
  n_sds = p_step(n_s, n_p, max_state);
  free(max_state);

  return n_sds;
}

int parity_from_p(unsigned int p, int n_s, int n_p, int* l_shell) {
/* Computes the parity eigenvalue +/- of a given SD p

  Input(s):
    unsigned int p: Slater determinant
    int n_s: number of single-particle states
    int n_p: number of particles
    int* l_shell: array of orbital angular momentum l values for each single-particle state

  Output(s):
    int pi_tot: total parity of the given SD
*/

  unsigned int q = n_choose_k(n_s, n_p) - p;
  int j_min = 1;
  int pi_tot = 1;
  for (int k = 0; k < n_p; k++) {
    for (int j = j_min; j <= n_s; j++) {
      if ((q >= n_choose_k(n_s - j, n_p - k)) && (q < n_choose_k(n_s - (j - 1), n_p - k))) {
        pi_tot *= pow(-1.0, l_shell[j - 1]);
        q -= n_choose_k(n_s - j, n_p  - k);
        j_min = j+1;
        break;
      }
    }
  }

  return pi_tot;
}

int w_from_p(unsigned int p, int n_s, int n_p, int* w_shell) {
  /* Returns the total weighted w of a given SD p-coefficient
  Uses the same algorithm as orbitals_from_p to reconstruct the occupied orbitals
  and simply sums the resulting w values

  Input(s): 
    unsigned int p: SD p-coefficient
    int n_s: number of single-particle states
    int n_p: number of particles
    int array w_shell: array of w values for each single-particle orbital
    
  Output(s): 
    int m_tot: total m_j value of the given SD
 */

  unsigned int q = n_choose_k(n_s, n_p) - p;
  int j_min = 1;
  int w_tot = 0;
  for (int k = 0; k < n_p; k++) {
    for (int j = j_min; j <= n_s; j++) {
      if ((q >= n_choose_k(n_s - j, n_p - k)) && (q < n_choose_k(n_s - (j - 1), n_p - k))) {
        w_tot += w_shell[j-1];
        q -= n_choose_k(n_s - j, n_p  - k);
        j_min = j+1;
        break;
      }
    }
  }
  
  return w_tot;
}


float m_from_p(unsigned int p, int n_s, int n_p, int* m_shell) {
  /* Returns the total magnetic angular momentum m_j of a given SD p-coefficient
  Uses the same algorithm as orbitals_from_p to reconstruct the occupied orbitals
  and simply sums the resulting m_j values

  Input(s): 
    unsigned int p: SD p-coefficient
    int n_s: number of single-particle states
    int n_p: number of particles
    int array m_shell: array of m_j values for each single-particle orbital
    
  Output(s): 
    int m_tot: total m_j value of the given SD
 */

  unsigned int q = n_choose_k(n_s, n_p) - p;
  int j_min = 1;
  float m_tot = 0;
  for (int k = 0; k < n_p; k++) {
    for (int j = j_min; j <= n_s; j++) {
      if ((q >= n_choose_k(n_s - j, n_p - k)) && (q < n_choose_k(n_s - (j - 1), n_p - k))) {
        m_tot += m_shell[j-1]/2.0;
        q -= n_choose_k(n_s - j, n_p  - k);
        j_min = j+1;
        break;
      }
    }
  }
  return m_tot;
}

float max_mj(int n_s, int n_p, int *m_shell) {
  /* Computes the maximum mj value of all SDs in the given basis
  
  Input(s):
    int n_s: number of single-particle states
    int n_p: number of particles
    int* m_shell: array of mj values for each single-particle state

  Outputs():
    int max_mj: maximum mj in the basis
*/
   
  int *found = (int*) calloc(n_s, sizeof(int));
  float max_j = 0;
  for (int i = 0; i < n_p; i++) {
    float cur_max_j = -100;
    int i_found;
    for (int j = 0; j < n_s; j++) {
      if (found[j]) {continue;}
      if (m_shell[j]/2.0 > cur_max_j) {
        cur_max_j = m_shell[j]/2.0;
        i_found = j;
      }
    }
    found[i_found] = 1;
    max_j += cur_max_j;
  }
  free(found);

  return max_j;
}

float min_w(int n_s, int n_p, int *w_shell) {
  /* Computes the minimum w value of all SDs in the given basis
  
  Input(s):
    int n_s: number of single-particle states
    int n_p: number of particles
    int* m_shell: array of mj values for each single-particle state

  Outputs():
    int min_mj: minimum mj in the basis
*/

  int *found = (int*) calloc(n_s, sizeof(int));
  int min_w = 0;
  for (int i = 0; i < n_p; i++) {
    int cur_min_w = 10000;
    int i_found;
    for (int j = 0; j < n_s; j++) {
      if (found[j]) {continue;}
      if (w_shell[j] < cur_min_w) {
        cur_min_w = w_shell[j];
        i_found = j;
      }
    }
    min_w += cur_min_w;
    found[i_found] = 1;
  }
  free(found);

  return min_w;
}


float min_mj(int n_s, int n_p, int *m_shell) {
  /* Computes the minimum mj value of all SDs in the given basis
  
  Input(s):
    int n_s: number of single-particle states
    int n_p: number of particles
    int* m_shell: array of mj values for each single-particle state

  Outputs():
    int min_mj: minimum mj in the basis
*/

  int *found = (int*) calloc(n_s, sizeof(int));
  float min_j = 0;
  for (int i = 0; i < n_p; i++) {
    float cur_min_j = 100;
    int i_found;
    for (int j = 0; j < n_s; j++) {
      if (found[j]) {continue;}
      if (m_shell[j]/2.0 < cur_min_j) {
        cur_min_j = m_shell[j]/2.0;
        i_found = j;
      }
    }
    min_j += cur_min_j;
    found[i_found] = 1;
  }
  free(found);

  return min_j;
}

int j_min_from_p(int n_s, int n_p, unsigned int p) {
/* Returns the label of the first occupied orbital
   Uses the standard Whitehead algorithm for reconstructing occupied states from p-coefficient

  Input(s):
    int n_s: number of single-particle states
    int n_p: number of particles
    unsigned int p: SD p-coefficient

  Output(s):
    int j_min: label of first occupied orbital

*/

  unsigned int q = n_choose_k(n_s, n_p) - p;
  
  for (int j = 1; j <= n_s; j++) {
    if ((q >= n_choose_k(n_s - j, n_p)) && (q < n_choose_k(n_s - (j - 1), n_p))) {
      return j;
    }
  }

  return 0;
}

unsigned long n_choose_k(unsigned int n, unsigned int k) {
/* Code to compute binomial coefficient
   Provides zero result when n < k (required for Whitehead algorithms)
*/
//  if (n < 0 || k < 0) {printf("Error %d %d\n", n, k);}
  unsigned int c = 0;
  if (k > n) {return c;}
  if (k == 0) {return 1;}
//  c = gsl_sf_choose(n, k);
//  return c;
  return ((n * n_choose_k(n-1, k - 1)) / k);
}

int test_n_choose_k() {
  int pass = 1;

  if (n_choose_k(0,0) != 1) {pass = 0;}
  if (n_choose_k(1,0) != 1) {pass = 0;}
  if (n_choose_k(0,1) != 0) {pass = 0;}
  if (n_choose_k(6,2) != 15) {pass = 0;}
  if (n_choose_k(25, 6) != 177100) {pass = 0;}
  if (n_choose_k(52, 27) != 477551179875952) {pass = 0;}

  if (pass == 1) {
    printf("Test n_choose_k: Pass\n");
  } else {
    printf("Test n_choose_k: Fail\n");
  }

  return pass;
}
