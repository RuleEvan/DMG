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
  unsigned int p = gsl_sf_choose(n_s, n_p);
  for (int i = 0; i < n_p; i++) {
    p -= n_choose_k(n_s-m_p[i], n_p - i);
  }
  return p;
}

void orbitals_from_p(unsigned int p, int n_s, int n_p, int* orbitals) {
/* Given a p-coefficient, reconstruct the orbitals occupied in
   the corresponding Slater determinant and store in orbitals array
   see Whitehead for details

  Input(s):
    int p: p-coefficient
    int n_s: number of single-particle states states
    int n_p: number of particles
    int array orbitals: array to be populated with orbital number
                        of each particle

  Output(s): None

*/
  unsigned int q = n_choose_k(n_s, n_p) - p;
  int j_min = 1;
  for (int k = 0; k < n_p; k++) {
    for (int j = j_min; j <= n_s; j++) {
      if ((q >= n_choose_k(n_s - j, n_p - k)) && (q < n_choose_k(n_s - (j - 1), n_p - k))) {
        orbitals[k] = j;
        q -= n_choose_k(n_s - j, n_p  - k);
        j_min = j+1;
        break;
      }
    }
  }
  return;
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
      // then we know that this state is unoccied and we can return zero
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

void get_m_pi_q(unsigned int p, int n_s, int n_p, int* n_shell, int* l_shell, int* jz_shell, float* m_j, int* parity, int* n_quanta, int j_min) {
  unsigned int q = n_choose_k(n_s, n_p) - p;
  *parity = 1;
  *m_j = 0;
  *n_quanta = 0;
  for (int k = 0; k < n_p; k++) {
    for (int j = j_min; j <= n_s; j++) {
      if ((q >= n_choose_k(n_s - j, n_p - k)) && (q < n_choose_k(n_s - (j - 1), n_p - k))) {
        *m_j += jz_shell[j - 1]/2.0;
        *parity *= pow(-1.0, l_shell[j - 1]);
        *n_quanta += 2*n_shell[j - 1] + l_shell[j - 1];
     
        q -= n_choose_k(n_s - j, n_p  - k);
        j_min = j+1;
        break;
      }
    }
  }

  return;
}

 
int get_num_quanta(unsigned int p, int n_s, int n_p, int* n_shell, int* l_shell) {
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
  int n_tot = 1;
  for (int k = 0; k < n_p; k++) {
    for (int j = j_min; j <= n_s; j++) {
      if ((q >= n_choose_k(n_s - j, n_p - k)) && (q < n_choose_k(n_s - (j - 1), n_p - k))) {
        n_tot += 2*n_shell[j - 1] + l_shell[j - 1];
        q -= n_choose_k(n_s - j, n_p  - k);
        j_min = j+1;
        break;
      }
    }
  }

  return n_tot;
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

int get_max_n_spec_q(int n_s, int n_p, int *n_shell, int* l_shell) {
  /* Computes the maximum mj value of all SDs in the given basis
  
  Input(s):
    int n_s: number of single-particle states
    int n_p: number of particles
    int* m_shell: array of mj values for each single-particle state

  Outputs():
    int max_mj: maximum mj in the basis
*/
   
  int *found = (int*) calloc(n_s, sizeof(int));
  int max_n_q = 0;
  for (int i = 0; i < n_p; i++) {
    int cur_max = 0;
    int i_found;
    for (int j = 0; j < n_s; j++) {
      if (found[j]) {continue;}
      int n_q = 2*n_shell[j] + l_shell[j];
      if (n_q > cur_max) {
        cur_max = n_q;
        i_found = j;
      }
    }
    found[i_found] = 1;
    max_n_q += cur_max;
  }
  free(found);

  return max_n_q;
}

int get_min_n_spec_q(int n_s, int n_p, int *n_shell, int* l_shell) {
  /* Computes the minimum mj value of all SDs in the given basis
  
  Input(s):
    int n_s: number of single-particle states
    int n_p: number of particles
    int* m_shell: array of mj values for each single-particle state

  Outputs():
    int min_mj: minimum mj in the basis
*/
  if (n_p - 2 <= 0) {return 0;}
  int *found = (int*) calloc(n_s, sizeof(int));
  int min_n_q = 0;
  for (int i = 0; i < n_p - 2; i++) {
    float cur_min = 1000;
    int i_found;
    for (int j = 0; j < n_s; j++) {
      if (found[j]) {continue;}
      int n_q = 2*n_shell[j] + l_shell[j];
      if (n_q < cur_min) {
        cur_min = n_q;
        i_found = j;
      }
    }
    min_n_q += cur_min;
    found[i_found] = 1;
  }
  free(found);

  return min_n_q;
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
/* Returns the label of the first occupited orbital
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

unsigned int n_choose_k(int n, int k) {
/* Wrapper to gsl code to compute binomial coefficient
   Provides zero result when n < k (required for Whitehead algorithms)
*/

  unsigned int c = 0;
  if (k > n) {return c;}
  c = gsl_sf_choose(n, k);

  return c;
}

/* **** Deprecated code ****

void orbitals_from_binary(int n_s, int n_p, unsigned int b, int* orbitals) {
  // Computes the orbital of each particle from the corresponding
  // binary integer
  int j = n_p - 1;
  for (int i = 0; i < n_s; i++) {
    if (b % (unsigned int) pow(2, i + 1) != 0) { 
      orbitals[j] = n_s - i;
      b -= pow(2, i);
      j--;
    }
  }
  return;
}

unsigned int bin_from_orbitals(int n_s, int n_p, int* orbitals) {
  unsigned int b = 0;
  for (int i = 0; i < n_p; i++) {
    b += pow(2, n_s - orbitals[i]);
  }
  return b;
}

unsigned int p_from_binary(int n_s, int n_p, unsigned int b) {
  // Computes the p-coefficient from the corresponding binary number
  int* orbitals = (int*) malloc(sizeof(int)*n_p);
  orbitals_from_binary(n_s, n_p, b, orbitals);
  unsigned int p = p_step(n_s, n_p, orbitals);
  free(orbitals);

  return p;
}


unsigned int bin_from_p(int n_s, int n_p, unsigned int p) {
  // Returns the binary SD from the corresponding p-value
  unsigned int q = n_choose_k(n_s, n_p) - p;
  int j_min = 1;
  unsigned int bin = 0;
  for (int k = 0; k < n_p; k++) {
    for (int j = j_min; j <= n_s; j++) {
      if ((q >= n_choose_k(n_s - j, n_p - k)) && (q < n_choose_k(n_s - (j - 1), n_p - k))) {
        bin += pow(2, n_s - j);
        q -= n_choose_k(n_s - j, n_p  - k);
        j_min = j+1;
        break;
      }
    }
  }
  return bin;
}

unsigned int bin_phase_from_p(int n_s, int n_p, unsigned int p, int n_op, int* phase) {
  // Computes the binary SD and corresponding phase from acting an
  // operator at position n_op
  unsigned int q = n_choose_k(n_s, n_p) - p;
  int j_min = 1;
  unsigned int bin = 0;
  *phase = 1;
  for (int k = 0; k < n_p; k++) {
    for (int j = j_min; j <= n_s; j++) {
      if ((q >= n_choose_k(n_s - j, n_p - k)) && (q < n_choose_k(n_s - (j - 1), n_p - k))) {
        bin += pow(2, n_s - j);
        q -= n_choose_k(n_s - j, n_p  - k);
        j_min = j+1;
        if (j < n_op) {*phase *= -1;}
        break;
      }
    }
  }
  return bin;
}


unsigned int a_dag_op(int n_s, int n_p, unsigned int p, int n_op, int *phase) {
  // Acts a creation operator on the orbital in position n_op
  *phase = 1;
  unsigned int bin = bin_phase_from_p(n_s, n_p, p, n_op, phase);
  unsigned int check = pow(2, n_s - n_op);
  if (bin & check) {return 0;}
  bin += check;
  unsigned int pp = p_from_binary(n_s, n_p + 1, bin);
  return pp;
}

unsigned int a_dag_a_op(int n_s, int n_p, unsigned int p, int n_a, int n_b, int *phase) {
  *phase = 1;
  unsigned int q = n_choose_k(n_s, n_p) - p;
  unsigned int pf = n_choose_k(n_s, n_p);
  int j_min = 1;
  //if (n_a > n_b) {*phase = -1;}
//  printf("n_a: %d n_b: %d\n", n_a, n_b);
  int a_found = 0;
  int b_found = 0;
  int k = 0;
  int kf = 0; // Tracks particles in new SD
  // k tracks particles in old SD
  for (int j = 1; j <= n_s; j++) {
    // Check if the state is present in the old SD
    if ((q >= n_choose_k(n_s - j, n_p - k)) && (q < n_choose_k(n_s - (j - 1), n_p - k))) {
      // Check if the state is to be added but not deleted first
      if ((j == n_a) && (j != n_b)) {return 0;}
      // Otherwise check if the state is to be deleted
      if (j != n_b) {
      // Found normal orbital
        if (kf >= n_p) {return 0;}
        // Adjust final SD using j and kf
        pf -= n_choose_k(n_s - j, n_p - kf);
        // Adjust old SD using j and k
        q -= n_choose_k(n_s - j, n_p  - k);
        if (j < n_a) {*phase *= -1;}
        if (j < n_b) {*phase *= -1;}
        // Increment kf
        kf++;
        k++;
      } else if (n_a == n_b) {
        // Found orbital where a=b
        a_found = 1;
        b_found = 1;
        if (kf >= n_p) {return 0;}
        q -= n_choose_k(n_s - j, n_p - k);
        pf -= n_choose_k(n_s - j, n_p - kf);
        kf++;
        k++;
      } else {
      //  Found orbital to delete
        b_found = 1;
        q -= n_choose_k(n_s - j, n_p - k);
        k++;
      }
    } else {
      // Unoccupied states reach this point
      // Check if state is to be deleted
      if (j == n_b) {return 0;}
      // Otherwise check if it is to be created
      if (j == n_a) {
      //  printf("Found orbital to add a: %d\n", j);
        a_found = 1;
        if (kf >= n_p) {return 0;}
        pf -= n_choose_k(n_s - j, n_p - kf);
        kf++;
      } 
    }
  }
  if (!(a_found && b_found)) {return 0;}
  if (kf != n_p) {return 0;}
  if (n_a == n_b) {pf = p;}

  return pf;
}

unsigned int a_op_b(int n_s, unsigned int b, int n_a, int *phase) {
  *phase = 1;
  unsigned int b_a = pow(2, n_s - n_a);
  if (!(b_a & b)) {return 0;}
  unsigned int bf = b - b_a;
  unsigned int bp = pow(2, n_s) - b_a;
  bp = bp & b;
  *phase = pow(-1, bit_count(bp));

  return bf;
}

unsigned int a2_op_b(int n_s, unsigned int b, int n_a, int n_b, int *phase) {
  *phase = 1;
  if (n_a == n_b) {return b;}
  unsigned int b_a = pow(2, n_s - n_a);
  unsigned int b_b = pow(2, n_s - n_b);
  if (b_a & b) {return 0;}
  unsigned int bf = b - b_b + b_a;

  unsigned int bp = 0;
  if (b_a > b_b) {
    bp = b_a - 2*b_b;
  } else {
    bp = b_b - 2*b_a;
  }
  bp = bp & b;
  *phase = (int) pow(-1.0, bit_count(bp));

  return bf;
}

unsigned int a4_op_b(int n_s, unsigned int b, int n_a, int n_b, int n_c, int n_d, int* phase) {
  *phase = 1;
  if ((n_a == n_c) && (n_b == n_d)) {return b;}
  if ((n_a == n_d) && (n_b == n_c)) {*phase = -1; return b;}

  unsigned int b_a = pow(2, n_s - n_a);
  unsigned int b_b = pow(2, n_s - n_b);
  unsigned int b_c = pow(2, n_s - n_c);
  unsigned int b_d = pow(2, n_s - n_d);

  unsigned int bf = b + b_a + b_b - b_c - b_d;
  unsigned int bp1, bp2;
  
  if (n_b == n_d) {
    if (b_a > b_c) {
      bp2 = b_a - 2*b_c;
    } else {
      bp2 = b_c - 2*b_a;
    }
    bp2 = bp2 & b;
    *phase = pow(-1, bit_count(bp2));
  } else {
    if (b_b == b_c) {
      bp1 = 0;
    } else if (b_c > b_b) {
      bp1 = b_c - 2*b_b;
    } else {
      bp1 = b_b - 2*b_c;
    }
    if (b_a == b_d) {
      bp2 = 0;
    } else if (b_d > b_a) {
      bp2 = b_d - 2*b_a;
    } else {
      bp2 = b_a - 2*b_d;
    }
    if (b_b != b_c) {bp1 = bp1 & b;}
    if (b_a != b_d) {bp2 = bp2 & (b + b_b - b_c);}
    *phase = pow(-1, 1 + bit_count(bp1 ^ bp2));
  }

  return bf;
}

int i_compare(const void * a, const void * b) {
  return (*(int*)a - *(int*)b);
}

int bit_count(unsigned int b) {
  int count = 0;
  while (b != 0) {
    b = b & (b - 1);
    count++;
  }
  return count;
}
unsigned int a4_op(int n_s, int n_p, unsigned int p, int n_a, int n_b, int n_c, int n_d, int *phase) {
  // Acts the operator a^(dag)(a)a^(dag)(b)a(d)a(c) on the given SD
  *phase = 1;
  if (n_c == n_d) {return 0;}
  if (n_a == n_b) {return 0;}
  unsigned int q = n_choose_k(n_s, n_p) - p;
  unsigned int pf = n_choose_k(n_s, n_p);
  // Determine phase factors
//  if (n_d > n_c) {*phase *= -1;}
//  if (n_b > n_c) {*phase *= -1;}
//  if (n_a > n_c) {*phase *= -1;}
//  if (n_b > n_d) {*phase *= -1;}
//  if (n_a > n_d) {*phase *= -1;}
  if (n_a > n_b && (n_b != n_d) && (n_a != n_d) && (n_a != n_c) && (n_b != n_c)) {*phase *= -1;}
  if (n_c > n_d && (n_b != n_d) && (n_a != n_d) && (n_a != n_c) && (n_b != n_c)) {*phase *= -1;} 
//  printf("n_a: %d n_b: %d\n", n_a, n_b);
  int a_found = 0;
  int b_found = 0;
  int c_found = 0;
  int d_found = 0;
  int kf = 0;
  int k = 0;
  if ((n_a == n_d) && (n_b == n_c)) {*phase *= -1;}
  for (int j = 1; j <= n_s; j++) {
   // Check if this state is occupied
   if ((q >= n_choose_k(n_s - j, n_p - k)) && (q < n_choose_k(n_s - (j - 1), n_p - k))) {
      // Check if either creation operator appears
      if ((j == n_a) || (j == n_b)) {
        // Check if the correponding annihilation operator appears
        if ((j != n_c) && (j != n_d)) {return 0;} // Attempting to fill occupied orbital
        //  Found orbital where a or b = c or d
        if (j == n_a) {a_found = 1;}
        else if (j == n_b) {b_found = 1;}
        if (j == n_c) {c_found = 1;}
        else if (j == n_d) {d_found = 1;}
        if (kf >= n_p) {return 0;} // Final state has too many particles
        if (n_a == n_c) {
          if (((n_b < n_a) && (n_d > n_a)) || ((n_b > n_a) && (n_d < n_a))) {*phase *= -1;}
        }
        if (n_a == n_d) {
          if (!(((n_b < n_a) && (n_c > n_a)) || ((n_b > n_a) && (n_c < n_a)))) {*phase *= -1;}
        }
        if (n_b == n_d) {
          if (((n_a < n_b) && (n_c > n_b)) || ((n_a > n_b) && (n_c < n_b))) {*phase *= -1;}
        }
        if (n_b == n_c) {
          if (!(((n_a < n_b) && (n_d > n_b)) || ((n_a > n_b) && (n_d < n_b)))) {*phase *= -1;}
        }
        q -= n_choose_k(n_s - j, n_p - k);
        pf -= n_choose_k(n_s - j, n_p - kf);
        k++;
        kf++;
      } else if ((j != n_c) && (j != n_d)) {
      // Found normal orbital
        if (kf >= n_p) {return 0;} // Final state has too many particles
        q -= n_choose_k(n_s - j, n_p  - k);
        pf -= n_choose_k(n_s - j, n_p - kf);
        if (j < n_a) {*phase *= -1;}
        if (j < n_b) {*phase *= -1;}
        if (j < n_c) {*phase *= -1;}
        if (j < n_d) {*phase *= -1;}
        kf++;
        k++;
      } else {
        //  Found orbital to delete
        if (j == n_c) {c_found = 1;}
        else if (j == n_d) {d_found = 1;}
        q -= n_choose_k(n_s - j, n_p - k);
        k++;
      }
    } else {
      // Below here orbitals are empty
      if ((j == n_c) || (j == n_d)) {return 0;}
      if ((j == n_a) || (j == n_b)) {
        // Found orbital to add
        if (j == n_a) {a_found = 1;}
        else if (j == n_b) {b_found = 1;}
        if (kf >= n_p) {return 0;} // Final state has too many particles
        pf -= n_choose_k(n_s - j, n_p - kf);
        kf++;
      } 
    }
  }
  if (!(a_found && b_found && c_found && d_found)) {return 0;}
  if ((kf != n_p) || (k != n_p)) {return 0;}
  return pf;
}



void generate_binomial_file() {
  FILE *out_file;
  out_file = fopen("fort.13", "w");
  for (int j = 0; j < 121; j++) {
    for (int i = 0; i < 121; i++) {
      fprintf(out_file, "%d\n", n_choose_k(i, j));
    }
  }
  fclose(out_file);
  
  return;
}

void initialize_orbitals(int* n_shell, int* l_shell, int* j_shell, int* m_shell) {
  int n_s = 0;
  for (int q = 0; q <= N_OSC_QUANTA; q++) {
    int n_max;
    if ((q % 2) == 0) {n_max = q/2;}
    else {n_max = (q - 1)/2;}
    for (int n = 0; n <= n_max; n++) {
      int l = q - 2*n;
      for (int jj = abs(2*l-1); jj <= 2*l + 1; jj += 2) {
        for (int mjj = -jj; mjj <= jj; mjj += 2) {
          n_shell[n_s] = n;
          l_shell[n_s] = l;
          j_shell[n_s] = jj;
          m_shell[n_s] = mjj; 
          n_s++;
        }
      }
    }
  }
  
  return;
}

void generate_single_particle_states() {
  unsigned state;
  state = 0;
  int n_s = (N_OSC_QUANTA + 1)*(N_OSC_QUANTA + 2)*(N_OSC_QUANTA + 3)/3; 
  printf("Number of single particle states: %d\n", n_s);
  int* n_shell = (int*)malloc(sizeof(int)*n_s);
  int* l_shell = (int*)malloc(sizeof(int)*n_s);
  int* j_shell = (int*)malloc(sizeof(int)*n_s);
  int* m_shell = (int*)malloc(sizeof(int)*n_s);
  initialize_orbitals(n_shell, l_shell, j_shell, m_shell);
  int* p_state = (int*)malloc(sizeof(int)*N_PROTON);
  for (int i = 0; i < N_PROTON; i++) {
    p_state[i] = n_s - (N_PROTON - i - 1);
  }
  int n_sds_p = p_step(n_s, N_PROTON, p_state);
  printf("proton SDs: %d\n", n_sds_p);
  int* n_state = (int*)malloc(sizeof(int)*N_NEUTRON);
  for (int i = 0; i < N_NEUTRON; i++) {
    n_state[i] = n_s - (N_NEUTRON - i - 1);
  }
  int n_sds_n = p_step(n_s, N_NEUTRON, n_state);
  printf("neutron SDs: %d\n", n_sds_n);
  int match = 0;
  for (int p1 = 0; p1 < n_sds_p; p1++) {
    int m_tot_p1 = m_from_p(p1, n_s, N_PROTON, m_shell);
    for (int n1 = 0; n1 < n_sds_n; n1++) {
      int m_tot_n1 = m_from_p(n1, n_s, N_NEUTRON, m_shell);
      if (m_tot_p1 + m_tot_n1 != M_BASIS) {continue;}
      match++;
    }
  }
  printf("%d\n", match);
  double *hamiltonian = (double*)malloc(sizeof(double)*match*match);
  for (int i = 0; i < match*match; i++) {hamiltonian[i] = 0.0;}
  int k,l;
  k = 0;
  for (int p1 = 0; p1 < n_sds_p; p1++) {
    int m_tot_p1 = m_from_p(p1, n_s, N_PROTON, m_shell);
    for (int n1 = 0; n1 < n_sds_n; n1++) {
      int m_tot_n1 = m_from_p(n1, n_s, N_NEUTRON, m_shell);
      if ((m_tot_p1 + m_tot_n1) != M_BASIS) {continue;}
      l = 0;
      for (int p2 = 0; p2 < n_sds_p; p2++) {
        int m_tot_p2 = m_from_p(p2, n_s, N_PROTON, m_shell);
        for (int n2 = 0; n2 < n_sds_n; n2++) {
          int m_tot_n2 = m_from_p(n2, n_s, N_NEUTRON, m_shell);
          if ((m_tot_p2 + m_tot_n2) != M_BASIS) {continue;}
          hamiltonian[k*match + l] = k*l;        
          l++;
        }
      }
      k++;
    }
  }
  lanczos_eigenvalues(hamiltonian, match, N_LANCZOS);  

  return;
}
*/
