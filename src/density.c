#include "density.h"

/* This file contains routines to generate one-body and two-body density matrices
   The input files are wavefunction files written by BIGSTICK
   The code relies on a factorization between proton and neutron Slater determinants
*/

extern int VERBOSE;
int FORMAT = I_FORMAT;

void two_body_density(speedParams *sp) {
  // Master routine for generating two-body density matrices
  // Read in data 
  char wfn_file_initial[100];
  strcpy(wfn_file_initial, sp->initial_file_base);
  strcat(wfn_file_initial, ".wfn");
  char wfn_file_final[100];
  strcpy(wfn_file_final, sp->final_file_base);
  strcat(wfn_file_final, ".wfn");
  char basis_file_initial[100];
  strcpy(basis_file_initial, sp->initial_file_base);
  strcat(basis_file_initial, ".bas");
  char basis_file_final[100];
  strcpy(basis_file_final, sp->final_file_base);
  strcat(basis_file_final, ".bas");
  wfnData *wd = read_binary_wfn_data(wfn_file_initial, wfn_file_final, basis_file_initial, basis_file_final);

  int j_op = sp->j_op;
  int t_op = sp->t_op;
  int ns = wd->n_shells;

  // Compute magnetic isospin and check that the calculation makes physical sense
  float mti = 0.5*(wd->n_proton_i - wd->n_neutron_i);
  float mtf = 0.5*(wd->n_proton_f - wd->n_neutron_f);
  float mt_op = mtf - mti;
  if (fabs(mt_op) > t_op) {printf("Error: operator iso-spin is insufficient to mediate a transition between the given nuclides.\n"); exit(0);}
  if (fabs(mt_op) > 2) {printf("Error: Change in magnetic isospin is too large to be mediated by two-nucleon operator.\n"); exit(0);}

  // Determine number of intermediate Slater determinants
  // after the action of one or two annihilation operators
  unsigned int n_sds_p_int1 = get_num_sds(wd->n_shells, wd->n_proton_f - 1);
  unsigned int n_sds_p_int2 = get_num_sds(wd->n_shells, wd->n_proton_f - 2);
  unsigned int n_sds_n_int1 = get_num_sds(wd->n_shells, wd->n_neutron_f - 1);
  unsigned int n_sds_n_int2 = get_num_sds(wd->n_shells, wd->n_neutron_f - 2);

  if (VERBOSE) {printf("Intermediate SDs: p1: %d n1: %d, p2: %d, n2: %d\n", n_sds_p_int1, n_sds_n_int1, n_sds_p_int2, n_sds_n_int2);}

  // Determine the min and max mj values of proton/neutron SDs
  float mj_min_p_i = min_mj(ns, wd->n_proton_i, wd->jz_shell);
  float mj_max_p_i = max_mj(ns, wd->n_proton_i, wd->jz_shell);
  float mj_min_n_i = min_mj(ns, wd->n_neutron_i, wd->jz_shell);
  float mj_max_n_i = max_mj(ns, wd->n_neutron_i, wd->jz_shell);

  // Account for the factor that m_p + m_n = 0 in determining min/max mj
  if (fabs(((int) (2*wd->j_nuc_i[0])) % 2) > pow(10, -3)) {
    // half-integer mj case
    mj_min_p_i = MAX(mj_min_p_i, -mj_max_n_i + 0.5);
    mj_max_p_i = MIN(mj_max_p_i, -mj_min_n_i + 0.5);
    mj_min_n_i = MAX(mj_min_n_i, -mj_max_p_i + 0.5);
    mj_max_n_i = MIN(mj_max_n_i, -mj_min_p_i + 0.5);
  } else {
    mj_min_p_i = MAX(mj_min_p_i, -mj_max_n_i);
    mj_max_p_i = MIN(mj_max_p_i, -mj_min_n_i);
    mj_min_n_i = MAX(mj_min_n_i, -mj_max_p_i);
    mj_max_n_i = MIN(mj_max_n_i, -mj_min_p_i);
  }

  // Determine number of sectors (mp, mn)
  int num_mj_p_i = mj_max_p_i - mj_min_p_i + 1;
  int num_mj_n_i = mj_max_n_i - mj_min_n_i + 1;


  float mj_min_p_f = min_mj(ns, wd->n_proton_f, wd->jz_shell);
  float mj_max_p_f = max_mj(ns, wd->n_proton_f, wd->jz_shell);
  float mj_min_n_f = min_mj(ns, wd->n_neutron_f, wd->jz_shell);
  float mj_max_n_f = max_mj(ns, wd->n_neutron_f, wd->jz_shell);

  // Account for the factor that m_p + m_n = 0 in determining min/max mj
  if (fabs(((int) (2*wd->j_nuc_i[0])) % 2) > pow(10, -3)) {
    // half-integer mj case
    mj_min_p_f = MAX(mj_min_p_f, -mj_max_n_f + 0.5);
    mj_max_p_f = MIN(mj_max_p_f, -mj_min_n_f + 0.5);
    mj_min_n_f = MAX(mj_min_n_f, -mj_max_p_f + 0.5);
    mj_max_n_f = MIN(mj_max_n_f, -mj_min_p_f + 0.5);
  } else { 
    mj_min_p_f = MAX(mj_min_p_f, -mj_max_n_f);
    mj_max_p_f = MIN(mj_max_p_f, -mj_min_n_f);
    mj_min_n_f = MAX(mj_min_n_f, -mj_max_p_f);
    mj_max_n_f = MIN(mj_max_n_f, -mj_min_p_f);
  }

  // Determine number of sectors (mp, mn)
  int num_mj_p_f = mj_max_p_f - mj_min_p_f + 1;
  int num_mj_n_f = mj_max_n_f - mj_min_n_f + 1;

if (VERBOSE) {
  printf("Number of initial proton sectors:  %d\n", num_mj_p_i);
  printf("Number of initial neutron sectors: %d\n", num_mj_n_i);

  printf("Number of final proton sectors:  %d\n", num_mj_p_f);
  printf("Number of final neutron sectors: %d\n", num_mj_n_f);

  printf("mj_min_p_i: %g\n", mj_min_p_i);
  printf("mj_max_p_i: %g\n", mj_max_p_i);
  printf("mj_min_p_f: %g\n", mj_min_p_i);
  printf("mj_max_p_f: %g\n", mj_max_p_f);
  printf("mj_min_n_i: %g\n", mj_min_n_i);
  printf("mj_max_n_i: %g\n", mj_max_n_i);
  printf("mj_min_n_f: %g\n", mj_min_n_f);
  printf("mj_max_n_f: %g\n", mj_max_n_f);

}

  // Allocate space for jump lists
  wf_list **p0_list_i = (wf_list**) calloc(2*num_mj_p_i, sizeof(wf_list*));
  wf_list **n0_list_i = (wf_list**) calloc(2*num_mj_n_i, sizeof(wf_list*));
  sd_list **p1_list_i = (sd_list**) calloc(2*ns*num_mj_p_i, sizeof(sd_list*));
  sd_list **n1_list_i = (sd_list**) calloc(2*ns*num_mj_n_i, sizeof(sd_list*));
  sd_list **p1_list_f = (sd_list**) calloc(2*ns*num_mj_p_i, sizeof(sd_list*));
  sd_list **n1_list_f = (sd_list**) calloc(2*ns*num_mj_n_i, sizeof(sd_list*));
  sd_list **p2_list_i = (sd_list**) calloc(2*ns*ns*num_mj_p_i, sizeof(sd_list*));
  sd_list **n2_list_i = (sd_list**) calloc(2*ns*ns*num_mj_n_i, sizeof(sd_list*));
  sd_list **p2_list_f = (sd_list**) calloc(2*ns*ns*num_mj_p_i, sizeof(sd_list*));
  sd_list **n2_list_f = (sd_list**) calloc(2*ns*ns*num_mj_n_i, sizeof(sd_list*));
  sd_list **p2_list_if = (sd_list**) calloc(2*ns*ns*num_mj_p_i, sizeof(sd_list*));
  sd_list **n2_list_if = (sd_list**) calloc(2*ns*ns*num_mj_n_i, sizeof(sd_list*));

  int* p1_array_f = (int*) calloc(ns*n_sds_p_int1, sizeof(int));
  int* p2_array_f = (int*) calloc(ns*ns*n_sds_p_int2, sizeof(int));
  int* n1_array_f = (int*) calloc(ns*n_sds_n_int1, sizeof(int));
  int* n2_array_f = (int*) calloc(ns*ns*n_sds_n_int2, sizeof(int));

  if (wd->w_max_i > 0 || wd->w_max_f > 0) {if (VERBOSE) {printf("Truncation scheme detected...generating truncated jumps...\n");}
    int w_min_p_i = min_w(ns, wd->n_proton_i, wd->w_shell);
    int w_min_n_i = min_w(ns, wd->n_neutron_i, wd->w_shell);
    int w_min_p_f = min_w(ns, wd->n_proton_f, wd->w_shell);
    int w_min_n_f = min_w(ns, wd->n_neutron_f, wd->w_shell);
    
    int w_max_p_i = wd->w_max_i - w_min_n_i;
    int w_max_n_i = wd->w_max_i - w_min_p_i;
    int w_max_p_f = wd->w_max_f - w_min_n_f;
    int w_max_n_f = wd->w_max_f - w_min_p_f;

    if (wd->same_basis) {
    // We assume that same basis implies that magnetic isospin is unchanged
      if (VERBOSE) {printf("Building jumps for delta_MT = 0\n");}
      build_two_body_jumps_dmt0_trunc(wd->n_shells, wd->n_proton_i, mj_min_p_i, mj_max_p_i, num_mj_p_i, wd->n_sds_p_i, n_sds_p_int1, n_sds_p_int2, p1_array_f, p2_array_f, p0_list_i, p1_list_i, p2_list_i, wd->jz_shell, wd->l_shell, wd->w_shell, w_max_p_i);
      build_two_body_jumps_dmt0_trunc(wd->n_shells, wd->n_neutron_i, mj_min_n_i, mj_max_n_i, num_mj_n_i, wd->n_sds_n_i, n_sds_n_int1, n_sds_n_int2, n1_array_f, n2_array_f, n0_list_i, n1_list_i, n2_list_i, wd->jz_shell, wd->l_shell, wd->w_shell, w_max_n_i);
    } else if (mt_op == 2) {
      if (VERBOSE) {printf("Building jumps for delta_MT = +2\n");}
      build_two_body_jumps_dmtp2_trunc(wd->n_shells, wd->n_proton_f, num_mj_p_i, mj_min_p_i, mj_max_p_i, num_mj_p_f, mj_min_p_f, mj_max_p_f, wd->n_sds_p_f, p2_list_f, wd->n_neutron_i, num_mj_n_i, mj_min_n_i, mj_max_n_i, num_mj_n_f, mj_min_n_f, mj_max_n_f, wd->n_sds_n_i, n2_list_i, wd->jz_shell, wd->l_shell, wd->w_shell, w_max_p_i, w_max_p_f, w_max_n_i, w_max_n_f);
    } else if (mt_op == -2) {
      if (VERBOSE) {printf("Building jumps for delta_Mt = -2\n");}
      build_two_body_jumps_dmtp2_trunc(wd->n_shells, wd->n_neutron_f, num_mj_n_i, mj_min_n_i, mj_max_n_i, num_mj_n_f, mj_min_n_f, mj_max_n_f, wd->n_sds_n_f, n2_list_f, wd->n_proton_i, num_mj_p_i, mj_min_p_i, mj_max_p_i, num_mj_p_f, mj_min_p_f, mj_max_p_f, wd->n_sds_p_i, p2_list_i, wd->jz_shell, wd->l_shell, wd->w_shell, w_max_n_i, w_max_n_f, w_max_p_i, w_max_p_f);
    } else if (mt_op == 1) {
      if (VERBOSE) {printf("Building jumps for delta_Mt = +1\n");}
      build_two_body_jumps_dmtp1_trunc(wd->n_shells, wd->n_proton_i, num_mj_p_i, mj_min_p_i, mj_max_p_i, wd->n_proton_f, num_mj_p_f, mj_min_p_f, mj_max_p_f, p1_list_i, p1_list_f, p2_array_f, wd->n_neutron_i, num_mj_n_i, mj_min_n_i, mj_max_n_i, wd->n_neutron_f, num_mj_n_f, mj_min_n_f, mj_max_n_f, n1_list_i, n2_list_i, n1_array_f, wd->jz_shell, wd->l_shell, wd->w_shell, w_max_p_i, w_max_p_f, w_max_n_i, w_max_n_f);    
    } else if (mt_op == -1) {
      if (VERBOSE) {printf("Building jumps for delta_Mt = -1\n");}
      build_two_body_jumps_dmtp1_trunc(wd->n_shells, wd->n_neutron_i, num_mj_n_i, mj_min_n_i, mj_max_n_i, wd->n_neutron_f, num_mj_n_f, mj_min_n_f, mj_max_n_f, n1_list_i, n1_list_f, n2_array_f, wd->n_proton_i, num_mj_p_i, mj_min_p_i, mj_max_p_i, wd->n_proton_f, num_mj_p_f, mj_min_p_f, mj_max_p_f, p1_list_i, p2_list_i, p1_array_f, wd->jz_shell, wd->l_shell, wd->w_shell, w_max_n_i, w_max_n_f, w_max_p_i, w_max_p_f);
    }
  } else {
    if (wd->same_basis) {
    // We assume that same basis implies that magnetic isospin is unchanged
      if (VERBOSE) {printf("Building jumps for delta_MT = 0\n");}
      build_two_body_jumps_dmt0_alt(wd->n_shells, wd->n_proton_i, mj_min_p_i, mj_max_p_i, num_mj_p_i, wd->n_sds_p_i, n_sds_p_int1, n_sds_p_int2, p1_array_f, p2_array_f, p0_list_i, p1_list_i, p2_list_i, p2_list_if, wd->jz_shell, wd->l_shell);
      build_two_body_jumps_dmt0_alt(wd->n_shells, wd->n_neutron_i, mj_min_n_i, mj_max_n_i, num_mj_n_i, wd->n_sds_n_i, n_sds_n_int1, n_sds_n_int2, n1_array_f, n2_array_f, n0_list_i, n1_list_i, n2_list_i, n2_list_if, wd->jz_shell, wd->l_shell);
    } else if (mt_op == 2) {
      if (VERBOSE) {printf("Building jumps for delta_MT = +2\n");}
      build_two_body_jumps_dmtp2(wd->n_shells, wd->n_proton_f, num_mj_p_i, mj_min_p_i, mj_max_p_i, num_mj_p_f, mj_min_p_f, mj_max_p_f, wd->n_sds_p_f, p2_list_f, wd->n_neutron_i, num_mj_n_i, mj_min_n_i, mj_max_n_i, num_mj_n_f, mj_min_n_f, mj_max_n_f, wd->n_sds_n_i, n2_list_i, wd->jz_shell, wd->l_shell);
    } else if (mt_op == -2) {
      if (VERBOSE) {printf("Building jumps for delta_Mt = -2\n");}
      build_two_body_jumps_dmtp2(wd->n_shells, wd->n_neutron_f, num_mj_n_i, mj_min_n_i, mj_max_n_i, num_mj_n_f, mj_min_n_f, mj_max_n_f, wd->n_sds_n_f, n2_list_f, wd->n_proton_i, num_mj_p_i, mj_min_p_i, mj_max_p_i, num_mj_p_f, mj_min_p_f, mj_max_p_f, wd->n_sds_p_i, p2_list_i, wd->jz_shell, wd->l_shell);
    } else if (mt_op == 1) {
      if (VERBOSE) {printf("Building jumps for delta_Mt = +1\n");}
      build_two_body_jumps_dmtp1(wd->n_shells, wd->n_proton_i, num_mj_p_i, mj_min_p_i, mj_max_p_i, wd->n_proton_f, num_mj_p_f, mj_min_p_f, mj_max_p_f, p1_list_i, p1_list_f, p2_array_f, wd->n_neutron_i, num_mj_n_i, mj_min_n_i, mj_max_n_i, wd->n_neutron_f, num_mj_n_f, mj_min_n_f, mj_max_n_f, n1_list_i, n2_list_i, n1_array_f, wd->jz_shell, wd->l_shell);    
    } else if (mt_op == -1) {
      if (VERBOSE) {printf("Building jumps for delta_Mt = -1\n");}
      build_two_body_jumps_dmtp1(wd->n_shells, wd->n_neutron_i, num_mj_n_i, mj_min_n_i, mj_max_n_i, wd->n_neutron_f, num_mj_n_f, mj_min_n_f, mj_max_n_f, n1_list_i, n1_list_f, n2_array_f, wd->n_proton_i, num_mj_p_i, mj_min_p_i, mj_max_p_i, wd->n_proton_f, num_mj_p_f, mj_min_p_f, mj_max_p_f, p1_list_i, p2_list_i, p1_array_f, wd->jz_shell, wd->l_shell);
    }
  }


  double* cg_fact = (double*) calloc(sp->n_trans, sizeof(double));
  FILE *out_file;
  char output_density_file[100];
  char output_log_file[100];
  char output_suffix[100];
  strcpy(output_log_file, sp->out_file_base);
  strcat(output_log_file, ".log");

  // Loop over transitions
  eigen_list* trans = sp->transition_list;
  int i_trans = 0;
  while (trans != NULL) {
    int psi_i = trans->eig_i;
    int psi_f = trans->eig_f;
    float ji = wd->j_nuc_i[psi_i];
    float ti = wd->t_nuc_i[psi_i];
    double cg_j = 0.0;
    double cg_t = 0.0;
    float jf = wd->j_nuc_f[psi_f];
    float tf = wd->t_nuc_f[psi_f];
    int ijf = 2*jf;
    float mj = 0.0;
    if (ijf % 2) {
      mj = 0.5;
    }
    cg_j = clebsch_gordan(j_op, ji, jf, 0, mj, mj);
    if (cg_j == 0.0) {trans = trans->next; continue;}
    cg_t = clebsch_gordan(t_op, ti, tf, mt_op, mti, mtf);
    if (cg_t == 0.0) {trans = trans->next; continue;}
    cg_j *= pow(-1.0, jf - ji)/sqrt(2*jf + 1);
    cg_t *= pow(-1.0, tf - ti)/sqrt(2*tf + 1);
    if (VERBOSE) {printf("Initial state: # %d J: %g T: %g Final state: # %d J: %g T: %g\n", psi_i + 1, ji, ti, psi_f + 1, jf, tf);}
    strcpy(output_density_file, sp->out_file_base);
    sprintf(output_suffix, "_J%d_T%d_%d_%d.dens", j_op, t_op, psi_i, psi_f);
    strcat(output_density_file, output_suffix);
    out_file = fopen(output_density_file, "w"); 
    cg_fact[i_trans] = cg_j*cg_t;
    i_trans++;
    trans = trans->next;
  }
  
  double* j_store = (double*) malloc(4*sizeof(double));
  double* density = (double*) calloc(sp->n_trans, sizeof(double));
  // Loop over orbital a
  for (int i_orb1 = 0; i_orb1 < wd->n_orbits; i_orb1++) {
    float j1 = wd->j_orb[i_orb1];
    // Loop over orbital b
    for (int i_orb2 = 0; i_orb2 <= i_orb1; i_orb2++) {
      float j2 = wd->j_orb[i_orb2];
      // Loop over orbital c
      for (int i_orb3 = 0; i_orb3 < wd->n_orbits; i_orb3++) {
        float j4 = wd->j_orb[i_orb3];
        // Loop over orbital d
        for (int i_orb4 = 0; i_orb4 <= i_orb3; i_orb4++) {
          //printf("%d %d %d %d\n", i_orb1, i_orb2, i_orb3, i_orb4);
          if (pow(-1.0, wd->l_orb[i_orb1] + wd->l_orb[i_orb2] + wd->l_orb[i_orb3] + wd->l_orb[i_orb4]) != 1.0) {continue;} // Parity
          float j3 = wd->j_orb[i_orb4];
          // Allocate storage for each j12 and j34
          int j_min_12 = abs(j1 - j2);
          int j_max_12 = j1 + j2;
          int j_dim_12 = j_max_12 - j_min_12 + 1;
          int j_min_34 = abs(j3 - j4);
          int j_max_34 = j3 + j4;
          int j_dim_34 = j_max_34 - j_min_34 + 1;
          int j_dim = j_dim_12*j_dim_34;
          if (j_op > j_max_12 + j_max_34) {continue;}
          int j_tot_min = j_max_12 + j_max_34;
          for (int j12 = j_min_12; j12 <= j_max_12; j12++) {
            for (int j34 = j_min_34; j34 <= j_max_34; j34++) {
              if (abs(j12 - j34) < j_tot_min) {j_tot_min = abs(j12 - j34);}
            }
          }
          if (j_op < j_tot_min) {continue;}
          j_store = realloc(j_store, sizeof(double)*4*j_dim*sp->n_trans);
          for (int k = 0; k < 4*j_dim*sp->n_trans; k++) {
            j_store[k] = 0.0;
          }

          // Loop over shells for orbit a
          for (int ia = 0; ia < 2*wd->n_shells; ia++) {
            float mt1 = 0.5;
            int a = ia;
            if (a >= ns) {a -= ns; mt1 -= 1;}
            if (wd->l_shell[a] != wd->l_orb[i_orb1]) {continue;}
            if (wd->n_shell[a] != wd->n_orb[i_orb1]) {continue;}
            if (wd->j_shell[a]/2.0 != j1) {continue;}
            float mj1 = wd->jz_shell[a]/2.0;
            // Loop over shells for orbit b
            for (int ib = 0; ib < 2*wd->n_shells; ib++) {
              if (ib == ia) {continue;}
              float mt2 = 0.5;
              int b = ib;
              if (b >= ns) {b -= ns; mt2 -= 1;}
              if (wd->l_shell[b] != wd->l_orb[i_orb2]) {continue;}
              if (wd->n_shell[b] != wd->n_orb[i_orb2]) {continue;}
              if (wd->j_shell[b]/2.0 != j2) {continue;}
              float mj2 = wd->jz_shell[b]/2.0;
              if ((i_orb1 == i_orb2) && (a < b) && (mt1 == mt2)) {continue;}
              for (int ic = 0; ic < 2*wd->n_shells; ic++) {
                float mt4 = 0.5;
                int c = ic;
                if (c >= ns) {c -= ns; mt4 -= 1;}
                if (wd->l_shell[c] != wd->l_orb[i_orb3]) {continue;}
                if (wd->n_shell[c] != wd->n_orb[i_orb3]) {continue;}
                if (wd->j_shell[c]/2.0 != j4) {continue;}
                float mj4 = wd->jz_shell[c]/2.0;
                // Loop over shells for orbit d
                for (int id = 0; id < 2*wd->n_shells; id++) {
                  if (ic == id) {continue;}
                  int d = id;
                  float mt3 = 0.5;
                  if (d >= ns) {d -= ns; mt3 -= 1;}
                  if (wd->l_shell[d] != wd->l_orb[i_orb4]) {continue;}
                  if (wd->n_shell[d] != wd->n_orb[i_orb4]) {continue;}
                  if (wd->j_shell[d]/2.0 != j3) {continue;}
                  float mj3 = wd->jz_shell[d]/2.0;
                  if ((i_orb3 == i_orb4) && (c < d) && (mt3 == mt4)) {continue;}
                  if (mj1 + mj2 != mj3 + mj4) {continue;}
                  if (mt1 + mt2 - mt3 - mt4 != mt_op) {continue;}
                  for (int i = 0; i < sp->n_trans; i++) {density[i] = 0.0;}
		  if (mt_op == 0.0) { // mt_op = 0
                    if ((mt1 == 0.5) && (mt2 == 0.5) && (mt3 == 0.5) && (mt4 == 0.5)) { // p^dag p^dag p p
                        trace_two_body_nodes_dmt0_40(a, b, c, d, num_mj_p_i, mj_min_p_i, num_mj_n_i, mj_min_n_i, n_sds_p_int2, p2_array_f, p2_list_i, n0_list_i, wd, 0, sp->transition_list, density);
                    } else if ((mt3 == -0.5) && (mt4 == -0.5) && (mt1 == -0.5) && (mt2 == -0.5)) { // n^dag n^dag n n
                        trace_two_body_nodes_dmt0_40(a, b, c, d, num_mj_n_i, mj_min_n_i, num_mj_p_i, mj_min_p_i, n_sds_n_int2, n2_array_f, n2_list_i, p0_list_i, wd, 1, sp->transition_list, density);
	            } else if ((mt1 == 0.5) && (mt2 == -0.5) && (mt3 == 0.5) && (mt4 == -0.5)) { // p^dag n^dag p n
                      trace_two_body_nodes_dmt0_22_alt(a, d, b, c, num_mj_p_i, mj_min_p_i, mj_min_p_f, mj_max_p_f, num_mj_n_i, mj_min_n_i, mj_min_n_f, mj_max_n_f, p2_list_if, n2_list_if, wd, sp->transition_list, density, wd->jz_shell);
                      for (int i = 0; i < sp->n_trans; i++) {density[i] *= -1.0;}
                    } else if ((mt1 == -0.5) && (mt2 == 0.5) && (mt3 == 0.5) && (mt4 == -0.5)) { // n^dag p^dag p n
                      trace_two_body_nodes_dmt0_22_alt(b, d, a, c, num_mj_p_i, mj_min_p_i, mj_min_p_f, mj_max_p_f, num_mj_n_i, mj_min_n_i, mj_min_n_f, mj_max_n_f, p2_list_if, n2_list_if, wd, sp->transition_list, density, wd->jz_shell);
                    } else if ((mt1 == 0.5) && (mt2 == -0.5) && (mt3 == -0.5) && (mt4 == 0.5)) { // p^dag n^dag n p
                      trace_two_body_nodes_dmt0_22_alt(a, c, b, d, num_mj_p_i, mj_min_p_i, mj_min_p_f, mj_max_p_f, num_mj_n_i, mj_min_n_i, mj_min_n_f, mj_max_n_f, p2_list_if, n2_list_if, wd, sp->transition_list, density, wd->jz_shell);
                    } else if ((mt1 == -0.5) && (mt2 == 0.5) && (mt3 == -0.5) && (mt4 == 0.5)) { // n^dag p^dag n p
                      trace_two_body_nodes_dmt0_22_alt(b, c, a, d, num_mj_p_i, mj_min_p_i, mj_min_p_f, mj_max_p_f, num_mj_n_i, mj_min_n_i, mj_min_n_f, mj_max_n_f, p2_list_if, n2_list_if, wd, sp->transition_list, density, wd->jz_shell);
                      for (int i = 0; i < sp->n_trans; i++) {density[i] *= -1.0;}
                    }
                   /* } else if ((mt1 == 0.5) && (mt2 == -0.5) && (mt3 == 0.5) && (mt4 == -0.5)) { // p^dag n^dag p n
                      trace_two_body_nodes_dmt0_22(a, d, b, c, num_mj_p_i, mj_min_p_i, mj_min_p_f, mj_max_p_f, num_mj_n_i, mj_min_n_i, mj_min_n_f, mj_max_n_f, n_sds_p_int1, n_sds_n_int1, p1_list_i, n1_list_i, p1_array_f, n1_array_f, wd, sp->transition_list, density, wd->jz_shell);
                      for (int i = 0; i < sp->n_trans; i++) {density[i] *= -1.0;}
                    } else if ((mt1 == -0.5) && (mt2 == 0.5) && (mt3 == 0.5) && (mt4 == -0.5)) { // n^dag p^dag p n
                      trace_two_body_nodes_dmt0_22(b, d, a, c, num_mj_p_i, mj_min_p_i, mj_min_p_f, mj_max_p_f, num_mj_n_i, mj_min_n_i, mj_min_n_f, mj_max_n_f, n_sds_p_int1, n_sds_n_int1, p1_list_i, n1_list_i, p1_array_f, n1_array_f, wd, sp->transition_list, density, wd->jz_shell);
                    } else if ((mt1 == 0.5) && (mt2 == -0.5) && (mt3 == -0.5) && (mt4 == 0.5)) { // p^dag n^dag n p
                      trace_two_body_nodes_dmt0_22(a, c, b, d, num_mj_p_i, mj_min_p_i, mj_min_p_f, mj_max_p_f, num_mj_n_i, mj_min_n_i, mj_min_n_f, mj_max_n_f, n_sds_p_int1, n_sds_n_int1,p1_list_i, n1_list_i, p1_array_f, n1_array_f, wd, sp->transition_list, density, wd->jz_shell);
                    } else if ((mt1 == -0.5) && (mt2 == 0.5) && (mt3 == -0.5) && (mt4 == 0.5)) { // n^dag p^dag n p
                      trace_two_body_nodes_dmt0_22(b, c, a, d, num_mj_p_i, mj_min_p_i, mj_min_p_f, mj_max_p_f, num_mj_n_i, mj_min_n_i, mj_min_n_f, mj_max_n_f, n_sds_p_int1, n_sds_n_int1, p1_list_i, n1_list_i, p1_array_f, n1_array_f, wd, sp->transition_list, density, wd->jz_shell);
                      for (int i = 0; i < sp->n_trans; i++) {density[i] *= -1.0;}
                    }*/

		  } else if (mt_op == 1.0) { // mt_op = +1
                    if ((mt1 == 0.5) && (mt2 == 0.5) && (mt3 == 0.5) && (mt4 == -0.5)) { // p^dag p^dag p n
                      trace_two_body_nodes_dmtp1_31(a, b, d, c, num_mj_p_i, mj_min_p_i, num_mj_n_i, mj_min_n_i, n_sds_p_int2, p2_array_f, p1_list_i, n1_list_i, wd, 0, sp->transition_list, density);
                      for (int i = 0; i < sp->n_trans; i++) {density[i] *= -1.0;}
                    } else if ((mt1 == 0.5) && (mt2 == 0.5) && (mt3 == -0.5) && (mt4 == 0.5)) { // p^dag p^dag n p
                      trace_two_body_nodes_dmtp1_31(a, b, c, d, num_mj_p_i, mj_min_p_i, num_mj_n_i, mj_min_n_i, n_sds_p_int2, p2_array_f, p1_list_i, n1_list_i, wd, 0, sp->transition_list, density);
                      for (int i = 0; i < sp->n_trans; i++) {density[i] *= 1.0;}
                    } else if ((mt1 == 0.5) && (mt2 == -0.5) && (mt3 == -0.5) && (mt4 == -0.5)) { // p^dag n^dag n n
                      trace_two_body_nodes_dmtp1_13(a, b, c, d, num_mj_n_i, mj_min_n_i, num_mj_p_i, mj_min_p_i, n_sds_n_int1, n1_array_f, n2_list_i, p1_list_f, wd, 1, sp->transition_list, density);
                      for (int i = 0; i < sp->n_trans; i++) {density[i] *= 1.0;}
                    } else if ((mt1 == -0.5) && (mt2 == 0.5) && (mt3 == -0.5) && (mt4 == -0.5)) { // n^dag p^dag n n
                      trace_two_body_nodes_dmtp1_13(b, a, c, d, num_mj_n_i, mj_min_n_i, num_mj_p_i, mj_min_p_i, n_sds_n_int1, n1_array_f, n2_list_i, p1_list_f, wd, 1, sp->transition_list, density);
                      for (int i = 0; i < sp->n_trans; i++) {density[i] *= -1.0;}
                    }
		  } else if (mt_op == -1) { // mt_op = -1
		    if ((mt1 == -0.5) && (mt2 == -0.5) && (mt3 == -0.5) && (mt4 == 0.5)) { // n^dag n^dag n p
                      trace_two_body_nodes_dmtp1_31(a, b, d, c, num_mj_n_i, mj_min_n_i, num_mj_p_i, mj_min_p_i, n_sds_n_int2, n2_array_f, n1_list_i, p1_list_i, wd, 1, sp->transition_list, density);
                      for (int i = 0; i < sp->n_trans; i++) {density[i] *= -1.0;}
                    } else if ((mt1 == -0.5) && (mt2 == -0.5) && (mt3 == 0.5) && (mt4 == -0.5)) { // n^dag n^dag p n
                      trace_two_body_nodes_dmtp1_31(a, b, c, d, num_mj_n_i, mj_min_n_i, num_mj_p_i, mj_min_p_i, n_sds_n_int2, n2_array_f, n1_list_i, p1_list_i, wd, 1, sp->transition_list, density);
                      for (int i = 0; i < sp->n_trans; i++) {density[i] *= 1.0;}
                    } else if ((mt1 == -0.5) && (mt2 == 0.5) && (mt3 == 0.5) && (mt4 == 0.5)) { // n^dag p^dag p p
                      trace_two_body_nodes_dmtp1_13(a, b, c, d, num_mj_p_i, mj_min_p_i, num_mj_n_i, mj_min_n_i, n_sds_p_int1, p1_array_f, p2_list_i, n1_list_f, wd, 0, sp->transition_list, density);
                      for (int i = 0; i < sp->n_trans; i++) {density[i] *= 1.0;}
                    } else if ((mt1 == 0.5) && (mt2 == -0.5) && (mt3 == 0.5) && (mt4 == 0.5)) { // p^dag n^dag p p
                      trace_two_body_nodes_dmtp1_13(b, a, c, d, num_mj_p_i, mj_min_p_i, num_mj_n_i, mj_min_n_i, n_sds_p_int1, p1_array_f, p2_list_i, n1_list_f, wd, 0, sp->transition_list, density);
                      for (int i = 0; i < sp->n_trans; i++) {density[i] *= -1.0;}
                    }
                  } else if (mt_op == 2) {  // p^dag p^dag n n
                      trace_two_body_nodes_dmtp2(a, b, d, c, num_mj_p_i, mj_min_p_i, num_mj_n_i, mj_min_n_i, n2_list_i, p2_list_f, wd, 0, sp->transition_list, density);
                  } else if (mt_op == -2) { // n^dag n^dag p p
                      trace_two_body_nodes_dmtp2(a, b, d, c, num_mj_n_i, mj_min_n_i, num_mj_p_i, mj_min_p_i, p2_list_i, n2_list_f, wd, 1, sp->transition_list, density);
                  } else {printf("Error in mt_op: %d\n", mt_op); exit(0);}
                  //printf("%g %g %g %g %g %g %g %g\n", mj1, mj2, mj3, mj4, mt1, mt2, mt3, mt4); 
                  for (int j12 = j_min_12; j12 <= j_max_12; j12++) {
                    if ((mj1 + mj2 > j12) || (mj1 + mj2 < -j12)) {continue;}
                    float cg_j12 = clebsch_gordan(j1, j2, j12, mj1, mj2, mj1 + mj2);
                    if (cg_j12 == 0.0) {continue;}
                    for (int t12 = 0; t12 <= 1; t12++) {
                      if ((mt1 + mt2 > t12) || (mt1 + mt2 < -t12)) {continue;}
                      float cg_t12 = clebsch_gordan(0.5, 0.5, t12, mt1, mt2, mt1 + mt2);
                      for (int j34 = j_min_34; j34 <= j_max_34; j34++) {
                        if ((mj3 + mj4 > j34) || (mj3 + mj4 < -j34)){continue;}
                        float cg_j34 = clebsch_gordan(j4, j3, j34, mj4, mj3, mj3 + mj4);
                        if (cg_j34 == 0.0) {continue;}
                        float cg_jop = clebsch_gordan(j_op, j34, j12, 0, mj3 + mj4, mj1 + mj2);
                        if (cg_jop == 0.0) {continue;}
                        for (int t34 = 0; t34 <= 1; t34++) {
                          if ((mt3 + mt4 > t34) || (mt3 + mt4 < -t34)) {continue;}
                          float cg_t34 = clebsch_gordan(0.5, 0.5, t34, mt4, mt3, mt3 + mt4);
                          if (cg_t34 == 0.0) {continue;}
                          float cg_top = clebsch_gordan(t_op, t34, t12, mt_op, mt3 + mt4, mt1 + mt2);
                          if (cg_top == 0.0) {continue;}
                          double d2 = cg_j12*cg_j34*cg_jop;
                          d2 *= cg_t12*cg_t34*cg_top;
                          d2 *= 0.25*pow(-1.0, -j34 + j12 - t34 + t12)/sqrt((2*j12 + 1)*(2*t12 + 1));
                          if (i_orb1 == i_orb2) {d2 *= sqrt(2);}
                          if (i_orb3 == i_orb4) {d2 *= sqrt(2);}
                          if (d2 == 0.0) {continue;}

                          for (int i = 0; i < sp->n_trans; i++) {
                            j_store[i + sp->n_trans*(t12 + 2*(t34 + 2*((j12 - j_min_12)*j_dim_34 + (j34 - j_min_34))))] += density[i]*d2/cg_fact[i];
                            if ((i_orb1 == i_orb2) && (mt1 == mt2)) {
                              j_store[i + sp->n_trans*(t12 + 2*(t34 + 2*((j12 - j_min_12)*j_dim_34 + (j34 - j_min_34))))] += pow(-1.0, j1 + j2 - j12 - t12)*density[i]*d2/cg_fact[i];
                            }
                            if ((i_orb3 == i_orb4) && (mt3 == mt4)) {
                              j_store[i + sp->n_trans*(t12 + 2*(t34 + 2*((j12 - j_min_12)*j_dim_34 + (j34 - j_min_34))))] += pow(-1.0, j3 + j4 - j34 - t34)*density[i]*d2/cg_fact[i];
                            }
                            if ((i_orb1 == i_orb2) && (mt1 == mt2) && (i_orb3 == i_orb4) && (mt3 == mt4)) {
                              j_store[i + sp->n_trans*(t12 + 2*(t34 + 2*((j12 - j_min_12)*j_dim_34 + (j34 - j_min_34))))] += pow(-1.0, j1 + j2 + j3 + j4 - j12 - j34 - t12 - t34)*density[i]*d2/cg_fact[i];
                            }
                          }
                        }
                      }
                    }
                  }
                }
              }
            }
          } 
          trans = sp->transition_list;
          i_trans = 0;
          while (trans != NULL) {
            int psi_i = trans->eig_i;
            int psi_f = trans->eig_f;
            strcpy(output_density_file, sp->out_file_base);
            sprintf(output_suffix, "_J%d_T%d_%d_%d.dens", j_op, t_op, psi_i, psi_f);
            strcat(output_density_file, output_suffix);
            out_file = fopen(output_density_file, "a"); 

            for (int t12 = 0; t12 <= 1; t12++) {
              for (int t34 = 0; t34 <= 1; t34++) {
                for (int ij12 = 0; ij12 < j_dim_12; ij12++) {
                  for (int ij34 = 0; ij34 < j_dim_34; ij34++) {
                    if (fabs(j_store[i_trans + sp->n_trans*(t12 + 2*(t34 + 2*(ij12*j_dim_34 + ij34)))]) < pow(10, -16)) {continue;}

                    fprintf(out_file, "%d,%g,%d,%g,%d,%d,%d,%g,%d,%g,%d,%d,%g\n", 2*wd->n_orb[i_orb1] + wd->l_orb[i_orb1], 2*wd->j_orb[i_orb1], 2*wd->n_orb[i_orb2] + wd->l_orb[i_orb2], 2*wd->j_orb[i_orb2], 2*(ij12 + j_min_12), 2*t12, 2*wd->n_orb[i_orb3] + wd->l_orb[i_orb3], 2*wd->j_orb[i_orb3], 2*wd->n_orb[i_orb4] + wd->l_orb[i_orb4], 2*wd->j_orb[i_orb4], 2*(ij34 + j_min_34), 2*t34, j_store[i_trans + sp->n_trans*(t12 + 2*(t34 + 2*(ij12*j_dim_34 + ij34)))]);

                    if (i_orb1 != i_orb2) { 
                      fprintf(out_file, "%d,%g,%d,%g,%d,%d,%d,%g,%d,%g,%d,%d,%g\n", 2*wd->n_orb[i_orb2] + wd->l_orb[i_orb2], 2*wd->j_orb[i_orb2], 2*wd->n_orb[i_orb1] + wd->l_orb[i_orb1], 2*wd->j_orb[i_orb1], 2*(ij12 + j_min_12), 2*t12, 2*wd->n_orb[i_orb3] + wd->l_orb[i_orb3], 2*wd->j_orb[i_orb3], 2*wd->n_orb[i_orb4] + wd->l_orb[i_orb4], 2*wd->j_orb[i_orb4], 2*(ij34 + j_min_34), 2*t34, pow(-1.0, j1 + j2 - ij12 - j_min_12 - t12)*j_store[i_trans + sp->n_trans*(t12 + 2*(t34 + 2*(ij12*j_dim_34 + ij34)))]);
                    }
                    if (i_orb3 != i_orb4) {
                        fprintf(out_file, "%d,%g,%d,%g,%d,%d,%d,%g,%d,%g,%d,%d,%g\n", 2*wd->n_orb[i_orb1] + wd->l_orb[i_orb1], 2*wd->j_orb[i_orb1], 2*wd->n_orb[i_orb2] + wd->l_orb[i_orb2], 2*wd->j_orb[i_orb2], 2*(ij12 + j_min_12), 2*t12, 2*wd->n_orb[i_orb4] + wd->l_orb[i_orb4], 2*wd->j_orb[i_orb4], 2*wd->n_orb[i_orb3] + wd->l_orb[i_orb3], 2*wd->j_orb[i_orb3], 2*(ij34 + j_min_34), 2*t34, pow(-1.0, j3 + j4 - ij34 - j_min_34 - t34)*j_store[i_trans + sp->n_trans*(t12 + 2*(t34 + 2*(ij12*j_dim_34 + ij34)))]);
                    }
                    if ((i_orb1 != i_orb2) && (i_orb3 != i_orb4)) {
                      fprintf(out_file, "%d,%g,%d,%g,%d,%d,%d,%g,%d,%g,%d,%d,%g\n", 2*wd->n_orb[i_orb2] + wd->l_orb[i_orb2], 2*wd->j_orb[i_orb2], 2*wd->n_orb[i_orb1] + wd->l_orb[i_orb1], 2*wd->j_orb[i_orb1], 2*(ij12 + j_min_12), 2*t12, 2*wd->n_orb[i_orb4] + wd->l_orb[i_orb4], 2*wd->j_orb[i_orb4], 2*wd->n_orb[i_orb3] + wd->l_orb[i_orb3], 2*wd->j_orb[i_orb3], 2*(ij34 + j_min_34), 2*t34, pow(-1.0, j1 + j2 + j3 + j4 - ij12 - j_min_12 - ij34 - j_min_34 - t12 - t34)*j_store[i_trans + sp->n_trans*(t12 + 2*(t34 + 2*(ij12*j_dim_34 + ij34)))]);

                    }
                  }
                }
              }
            }
            trans = trans->next;
            i_trans++;
            fclose(out_file);
          }
        }
      }
    }
  } 
  free(j_store); 

  return;
}


void one_body_density(speedParams* sp) {
  // Reads in BIGSTICK basis/wavefunction (.trwfn) files along with
  // orbit definition (.sp) files and constructs the one-body density matrices
  // for each initial and final state eigenfunction
  // Initial and final wave functions must share the same orbit file
  // J_op and T_op are the total spin and isospin of the operator

  // Read in data 
  char wfn_file_initial[100];
  strcpy(wfn_file_initial, sp->initial_file_base);
  strcat(wfn_file_initial, ".wfn");
  char wfn_file_final[100];
  strcpy(wfn_file_final, sp->final_file_base);
  strcat(wfn_file_final, ".wfn");
  char basis_file_initial[100];
  strcpy(basis_file_initial, sp->initial_file_base);
  strcat(basis_file_initial, ".bas");
  char basis_file_final[100];
  strcpy(basis_file_final, sp->final_file_base);
  strcat(basis_file_final, ".bas");
  wfnData *wd = read_binary_wfn_data(wfn_file_initial, wfn_file_final, basis_file_initial, basis_file_final);

  int ns = wd->n_shells;

  // Determine the number of intermediate proton SDs
  int n_sds_p_int = get_num_sds(wd->n_shells, wd->n_proton_f - 1);
  int n_sds_n_int = get_num_sds(wd->n_shells, wd->n_neutron_f - 1);
  
  // Determine the min and max mj values of proton/neutron SDs
  float mj_min_p_i = min_mj(ns, wd->n_proton_i, wd->jz_shell);
  float mj_max_p_i = max_mj(ns, wd->n_proton_i, wd->jz_shell);
  float mj_min_n_i = min_mj(ns, wd->n_neutron_i, wd->jz_shell);
  float mj_max_n_i = max_mj(ns, wd->n_neutron_i, wd->jz_shell);

  // Account for the factor that m_p + m_n = 0 in determining min/max mj
  if (fabs(((int) (2*wd->j_nuc_i[0])) % 2) > pow(10, -3)) {
    // half-integer mj case
    mj_min_p_i = MAX(mj_min_p_i, -mj_max_n_i + 0.5);
    mj_max_p_i = MIN(mj_max_p_i, -mj_min_n_i + 0.5);
    mj_min_n_i = MAX(mj_min_n_i, -mj_max_p_i + 0.5);
    mj_max_n_i = MIN(mj_max_n_i, -mj_min_p_i + 0.5);
  } else {
    mj_min_p_i = MAX(mj_min_p_i, -mj_max_n_i);
    mj_max_p_i = MIN(mj_max_p_i, -mj_min_n_i);
    mj_min_n_i = MAX(mj_min_n_i, -mj_max_p_i);
    mj_max_n_i = MIN(mj_max_n_i, -mj_min_p_i);
  }

  // Determine number of sectors (mp, mn)
  int num_mj_p_i = mj_max_p_i - mj_min_p_i + 1;
  int num_mj_n_i = mj_max_n_i - mj_min_n_i + 1;

  float mj_min_p_f = min_mj(ns, wd->n_proton_f, wd->jz_shell);
  float mj_max_p_f = max_mj(ns, wd->n_proton_f, wd->jz_shell);
  float mj_min_n_f = min_mj(ns, wd->n_neutron_f, wd->jz_shell);
  float mj_max_n_f = max_mj(ns, wd->n_neutron_f, wd->jz_shell);

  // Account for the factor that m_p + m_n = 0 in determining min/max mj
  if (fabs(((int) (2*wd->j_nuc_i[0])) % 2) > pow(10, -3)) {
    // half-integer mj case
    mj_min_p_f = MAX(mj_min_p_f, -mj_max_n_f + 0.5);
    mj_max_p_f = MIN(mj_max_p_f, -mj_min_n_f + 0.5);
    mj_min_n_f = MAX(mj_min_n_f, -mj_max_p_f + 0.5);
    mj_max_n_f = MIN(mj_max_n_f, -mj_min_p_f + 0.5);
  } else { 
    mj_min_p_f = MAX(mj_min_p_f, -mj_max_n_f);
    mj_max_p_f = MIN(mj_max_p_f, -mj_min_n_f);
    mj_min_n_f = MAX(mj_min_n_f, -mj_max_p_f);
    mj_max_n_f = MIN(mj_max_n_f, -mj_min_p_f);
  }

  // Determine number of sectors (mp, mn)
  int num_mj_p_f = mj_max_p_f - mj_min_p_f + 1;
  int num_mj_n_f = mj_max_n_f - mj_min_n_f + 1;

  float mti = 0.5*(wd->n_proton_i - wd->n_neutron_i);
  float mtf = 0.5*(wd->n_proton_f - wd->n_neutron_f);
  float mt_op = mtf - mti;
  if (fabs(mt_op > 1)) {printf("Error: change in magnetic isospin is too large to be mediated by single-nucleon operator.\n"); exit(0);}

  // Allocate space for jump lists
  wf_list **p0_list_i = (wf_list**) calloc(2*num_mj_p_i, sizeof(wf_list*));
  wf_list **n0_list_i = (wf_list**) calloc(2*num_mj_n_i, sizeof(wf_list*));
  sd_list **p1_list_i = (sd_list**) calloc(2*ns*num_mj_p_i, sizeof(sd_list*));
  sd_list **n1_list_i = (sd_list**) calloc(2*ns*num_mj_n_i, sizeof(sd_list*));
  sd_list **p1_list_f = (sd_list**) calloc(2*ns*num_mj_p_i, sizeof(sd_list*));
  sd_list **n1_list_f = (sd_list**) calloc(2*ns*num_mj_n_i, sizeof(sd_list*));
  int* p1_array_f = (int*) calloc(ns*n_sds_p_int, sizeof(int));
  int* n1_array_f = (int*) calloc(ns*n_sds_n_int, sizeof(int));
  
  if (wd->w_max_i > 0 || wd->w_max_f > 0) {if (VERBOSE) {printf("Truncation scheme detected...generating truncated jumps...\n");}
    int w_min_p_i = min_w(ns, wd->n_proton_i, wd->w_shell);
    int w_min_n_i = min_w(ns, wd->n_neutron_i, wd->w_shell);
    int w_min_p_f = min_w(ns, wd->n_proton_f, wd->w_shell);
    int w_min_n_f = min_w(ns, wd->n_neutron_f, wd->w_shell);
    
    int w_max_p_i = wd->w_max_i - w_min_n_i;
    int w_max_n_i = wd->w_max_i - w_min_p_i;
    int w_max_p_f = wd->w_max_f - w_min_n_f;
    int w_max_n_f = wd->w_max_f - w_min_p_f;


    if (wd->same_basis) {
      if (VERBOSE) {printf("Building initial and final state proton jumps...\n");}
      build_one_body_jumps_dmt0_trunc(wd->n_shells, wd->n_proton_i, mj_min_p_i, mj_max_p_i, num_mj_p_i, wd->n_sds_p_i, n_sds_p_int, p1_array_f, p0_list_i, p1_list_i, wd->jz_shell, wd->l_shell, wd->w_shell, w_max_p_i); 

      if (VERBOSE) {printf("Building initial and final state neutron jumps...\n");}
      build_one_body_jumps_dmt0_trunc(wd->n_shells, wd->n_neutron_i, mj_min_n_i, mj_max_n_i, num_mj_n_i, wd->n_sds_n_i, n_sds_n_int, n1_array_f, n0_list_i, n1_list_i, wd->jz_shell, wd->l_shell, wd->w_shell, w_max_n_i); 
    } else if (mt_op == 1) {
      if (VERBOSE) {printf("Delta MT = +1\n");}
      build_one_body_jumps_dmtp1_trunc(wd->n_shells, wd->n_proton_f, num_mj_p_i, mj_min_p_i, mj_max_p_i, num_mj_p_f, mj_min_p_f, mj_max_p_f, p1_list_f, wd->n_neutron_i, num_mj_n_i, mj_min_n_i, mj_max_n_i, n1_list_i, wd->jz_shell, wd->l_shell, wd->w_shell, w_max_p_i, w_max_p_f, w_max_n_i, w_max_n_f);
    } else if (mt_op == -1) {
      if (VERBOSE) {printf("Delta MT = -1\n");}
      build_one_body_jumps_dmtp1_trunc(wd->n_shells, wd->n_neutron_f, num_mj_n_i, mj_min_n_i, mj_max_n_i, num_mj_n_f, mj_min_n_f, mj_max_n_f, n1_list_f, wd->n_proton_i, num_mj_p_i, mj_min_p_i, mj_max_p_i, p1_list_i, wd->jz_shell, wd->l_shell, wd->w_shell, w_max_n_i, w_max_n_f, w_max_p_i, w_max_p_f);
    } else {printf("Error in delta MT\n"); exit(0);} 
  } else {
  if (wd->same_basis) {
      if (VERBOSE) {printf("Building initial and final state proton jumps...\n");}
      build_one_body_jumps_dmt0(wd->n_shells, wd->n_proton_i, mj_min_p_i, mj_max_p_i, num_mj_p_i, wd->n_sds_p_i, n_sds_p_int, p1_array_f, p0_list_i, p1_list_i, wd->jz_shell, wd->l_shell); 

      if (VERBOSE) {printf("Building initial and final state neutron jumps...\n");}
      build_one_body_jumps_dmt0(wd->n_shells, wd->n_neutron_i, mj_min_n_i, mj_max_n_i, num_mj_n_i, wd->n_sds_n_i, n_sds_n_int, n1_array_f, n0_list_i, n1_list_i, wd->jz_shell, wd->l_shell); 
    } else if (mt_op == 1) {
      if (VERBOSE) {printf("Delta MT = +1\n");}
      build_one_body_jumps_dmtp1(wd->n_shells, wd->n_proton_f, num_mj_p_i, mj_min_p_i, mj_max_p_i, num_mj_p_f, mj_min_p_f, mj_max_p_f, p1_list_f, wd->n_neutron_i, num_mj_n_i, mj_min_n_i, mj_max_n_i, n1_list_i, wd->jz_shell, wd->l_shell);
    } else if (mt_op == -1) {
      if (VERBOSE) {printf("Delta MT = -1\n");}
      build_one_body_jumps_dmtp1(wd->n_shells, wd->n_neutron_f, num_mj_n_i, mj_min_n_i, mj_max_n_i, num_mj_n_f, mj_min_n_f, mj_max_n_f, n1_list_f, wd->n_proton_i, num_mj_p_i, mj_min_p_i, mj_max_p_i, p1_list_i, wd->jz_shell, wd->l_shell);
    } else {printf("Error in delta MT\n"); exit(0);} 
  }
  // Loop over transitions
  eigen_list* trans = sp->transition_list;
  eigen_list *new_trans = malloc(sizeof(*new_trans));
  new_trans = NULL;
  int i_trans = 0;
  while (trans != NULL) {
    int psi_i = trans->eig_i;
    int psi_f = trans->eig_f;
    float ji = wd->j_nuc_i[psi_i];
    float ti = wd->t_nuc_i[psi_i];
    double cg_j = 0.0;
    double cg_t = 0.0;
    float jf = wd->j_nuc_f[psi_f];
    float tf = wd->t_nuc_f[psi_f];
    double mji = 0;
    double mjf = 0;
    int ijf = (int) 2*jf;
    if (ijf % 2) {mjf = 0.5; mji = 0.5;}
    for (int j_op = (int) (fabs(jf - ji)); j_op <= (int) (ji + jf); j_op++) { 
      for (int t_op = MAX(0,abs(ti-tf)); t_op <= MIN(1,abs(ti+tf)); t_op++) {
	if (fabs(mt_op) > t_op) {continue;}
        cg_j = clebsch_gordan(j_op, ji, jf, 0, mji,mjf);
        if (fabs(cg_j) < pow(10, -6)) {if (VERBOSE) {printf("Warning: CG coefficient is zero for allowed transition: Ji = %g Jf = %g Jop = %d \nComputing the density matrix for this operator requires the raising/lowering operators (not currently supported.)\n This entry will be omitted even though there could be non-zero contributions.\n", ji, jf, j_op);} continue;}
        cg_t = clebsch_gordan(t_op, ti, tf, mt_op, mti, mtf);
        if (fabs(cg_t) < pow(10, -6)) {if (VERBOSE) {printf("Warning: CG coefficient is zero for allowed transition: Ti = %g Tf = %g Top = %d \nComputing the density matrix for this operator requires the raising/lowering operators (not currently supported.) \n This entry will be omitted even though there could be non-zero contributions.\n", ti, tf, t_op);} continue;}
        if (new_trans == NULL) {
          new_trans = create_eigen_node(psi_i, psi_f, NULL);
        } else {
          eigen_append(new_trans, psi_i, psi_f);
        }
        i_trans++;
      }
    }
    trans = trans->next;
  }
  double* cg_fact = (double*) calloc(i_trans, sizeof(double));
  int* j_op_arr = (int*) calloc(i_trans, sizeof(int));
  int* t_op_arr = (int*) calloc(i_trans, sizeof(int));
  int* ji_arr = (int*) calloc(i_trans, sizeof(int));
  int* ti_arr = (int*) calloc(i_trans, sizeof(int));
  int* jf_arr = (int*) calloc(i_trans, sizeof(int));
  int* tf_arr = (int*) calloc(i_trans, sizeof(int));
  double* ef_arr = (double*) calloc(i_trans, sizeof(double));
  double* ei_arr = (double*) calloc(i_trans, sizeof(double));

  trans = sp->transition_list;
  i_trans = 0;
  while (trans != NULL) {
    int psi_i = trans->eig_i;
    int psi_f = trans->eig_f;
    float ji = wd->j_nuc_i[psi_i];
    float ti = wd->t_nuc_i[psi_i];
    double cg_j = 0.0;
    double cg_t = 0.0;
    float jf = wd->j_nuc_f[psi_f];
    float tf = wd->t_nuc_f[psi_f];
    double mji = 0;
    double mjf = 0;
    int ijf = (int) 2*jf;
    if (ijf % 2) {mjf = 0.5; mji = 0.5;}
    for (int j_op = (int) abs(jf - ji); j_op <= (int) (ji + jf); j_op++) {
      for (int t_op = MAX(0,abs(ti-tf)); t_op <= MIN(1,ti+tf); t_op++) {
	if (fabs(mt_op) > t_op) {continue;}
        cg_j = clebsch_gordan(j_op, ji, jf, 0, mji,mjf);
        if (fabs(cg_j) < pow(10, -6)) {continue;}
        cg_t = clebsch_gordan(t_op, ti, tf, mt_op, mti, mtf);
        if (fabs(cg_t) < pow(10, -6)) {continue;}
        cg_j *= pow(-1.0, j_op - ji + jf)*sqrt(2*j_op + 1)/sqrt(2*jf + 1);
        cg_t *= pow(-1.0, t_op - ti + tf)*sqrt(2*t_op + 1)/sqrt(2*tf + 1);
        cg_fact[i_trans] = cg_j*cg_t;
        j_op_arr[i_trans] = j_op;
        t_op_arr[i_trans] = t_op;
        ji_arr[i_trans] = (int) (2*ji);
        jf_arr[i_trans] = (int) (2*jf);
        ti_arr[i_trans] = (int) (2*ti);
        tf_arr[i_trans] = (int) (2*tf);
        ei_arr[i_trans] = wd->e_nuc_i[psi_i];
        ef_arr[i_trans] = wd->e_nuc_f[psi_f];
        i_trans++;
      }
    }
    trans = trans->next;
  } 


  sp->transition_list = new_trans;
  sp->n_trans = i_trans;


  // Loop over final state orbits
  double *total = (double*) calloc(sp->n_trans*pow(wd->n_orbits, 2), sizeof(double));
  double *density = (double*) calloc(sp->n_trans, sizeof(double));
  for (int i_orb1 = 0; i_orb1 < wd->n_orbits; i_orb1++) {
    float j1 = wd->j_orb[i_orb1];
    // Loop over initial state orbits
    for (int i_orb2 = 0; i_orb2 < wd->n_orbits; i_orb2++) {
      float j2 = wd->j_orb[i_orb2];
      // Loop over initial state SDs
      for (int ia = 0; ia < 2*wd->n_shells; ia++) {
        float mt1 = 0.5;
        int a = ia;
        if (a >= ns) {a -= ns; mt1 -= 1;}
        if (wd->l_shell[a] != wd->l_orb[i_orb1]) {continue;}
        if (wd->n_shell[a] != wd->n_orb[i_orb1]) {continue;}
        if (wd->j_shell[a]/2.0 != j1) {continue;}
        float mj1 = wd->jz_shell[a]/2.0;
        // Loop over initial state shells
        for (int ib = 0; ib < 2*wd->n_shells; ib++) {
          float mt2 = 0.5;
          int b = ib;
          if (b >= ns) {b -= ns; mt2 -= 1;}
          if (wd->l_shell[b] != wd->l_orb[i_orb2]) {continue;}
          if (wd->n_shell[b] != wd->n_orb[i_orb2]) {continue;}
          if (wd->j_shell[b]/2.0 != j2) {continue;}
          float mj2 = wd->jz_shell[b]/2.0;
          for (int i = 0; i < sp->n_trans; i++) {
            density[i] = 0.0;
          }
          
	  if ((mt1 == 0.5) && (mt2 == 0.5)) {
            trace_one_body_nodes_dmt0(a, b, num_mj_p_i, mj_min_p_i, num_mj_n_i, mj_min_n_i, n_sds_p_int, p1_array_f, p1_list_i, n0_list_i, wd, 0, sp->transition_list, density);
          } else if ((mt1 == -0.5) && (mt2 == -0.5)) {
            trace_one_body_nodes_dmt0(a, b, num_mj_n_i, mj_min_n_i, num_mj_p_i, mj_min_p_i, n_sds_n_int, n1_array_f, n1_list_i, p0_list_i, wd, 1, sp->transition_list, density);
          } else if ((mt1 == 0.5) && (mt2 == -0.5)) {
            trace_one_body_nodes_dmtp1(a, b, num_mj_p_i, mj_min_p_i, num_mj_n_i, mj_min_n_i, n1_list_i, p1_list_f, wd, 0, sp->transition_list, density);
          } else {
            trace_one_body_nodes_dmtp1(a, b, num_mj_n_i, mj_min_n_i, num_mj_p_i, mj_min_p_i, p1_list_i, n1_list_f, wd, 1, sp->transition_list, density);
          }
          for (int i = 0; i < sp->n_trans; i++) {
            int j_op = j_op_arr[i];
	    int t_op = t_op_arr[i];
            float d2 = clebsch_gordan(j1, j2, j_op, mj1, -mj2, 0);
            d2 *= clebsch_gordan(0.5, 0.5, t_op, mt1, -mt2, mt_op);
            d2 *= pow(-1.0, j2 - mj2 + 0.5 - mt2);
            if (d2 == 0.0) {continue;}

            total[i + sp->n_trans*(i_orb1 + i_orb2*wd->n_orbits)] += density[i]*d2/cg_fact[i];
	  }

        }
      }
    }
  } 
  FILE *out_file;
  char output_density_file[100];
  char output_log_file[100];
  char output_suffix[100];
  strcpy(output_log_file, sp->out_file_base);
  strcat(output_log_file, ".log");
  trans = sp->transition_list;
  i_trans = 0;
  int psi_i_prev = -1;
  int psi_f_prev = -1;
  if (FORMAT == 0) {
    while (trans != NULL) {
      int psi_i = trans->eig_i;
      int psi_f = trans->eig_f;
      int j_op = j_op_arr[i_trans];
      int t_op = t_op_arr[i_trans];
      int ji = ji_arr[i_trans];
      int jf = jf_arr[i_trans];
      int ti = ti_arr[i_trans];
      int tf = tf_arr[i_trans];
      if ((psi_i_prev == -1 && psi_f_prev == -1)) {
        strcpy(output_density_file, sp->out_file_base);
        sprintf(output_suffix, "_%d_%d.dens", psi_i, psi_f);
        strcat(output_density_file, output_suffix);
        out_file = fopen(output_density_file, "w+");
      } else if (psi_i != psi_i_prev || psi_f != psi_f_prev) {
        fclose(out_file);
        strcpy(output_density_file, sp->out_file_base);
        sprintf(output_suffix, "_%d_%d.dens", psi_i, psi_f);
        strcat(output_density_file, output_suffix);
        out_file = fopen(output_density_file, "w+");
      }
      fprintf(out_file, "    %d   %d   %d   %d   %d   %d\n", jf, tf, ji, ti, 2*j_op, 2*t_op);  
      for (int i_orb1 = 0; i_orb1 < wd->n_orbits; i_orb1++) {
        for (int i_orb2 = 0; i_orb2 < wd->n_orbits; i_orb2++) {
          if (fabs(total[i_trans + sp->n_trans*(i_orb1 + i_orb2*wd->n_orbits)]) > pow(10, -12)) {
            fprintf(out_file, "        %d       %g       %d        %g       %g\n", 2*wd->n_orb[i_orb1] + wd->l_orb[i_orb1], 2*wd->j_orb[i_orb1], 2*wd->n_orb[i_orb2] + wd->l_orb[i_orb2], 2*wd->j_orb[i_orb2], total[i_trans + sp->n_trans*(i_orb1 + i_orb2*wd->n_orbits)]);
          }
        }
      }
      fprintf(out_file, "       -1\n");
      i_trans++;
      psi_i_prev = psi_i;
      psi_f_prev = psi_f;
      trans = trans->next;
    }    
    fclose(out_file);
  } else if (FORMAT == 1) {

    while (trans != NULL) {
      int psi_i = trans->eig_i;
      int psi_f = trans->eig_f;
      int j_op = j_op_arr[i_trans];
      int t_op = t_op_arr[i_trans];
      int ji = ji_arr[i_trans];
      int jf = jf_arr[i_trans];
      int ti = ti_arr[i_trans];
      int tf = tf_arr[i_trans];
      double ei = ei_arr[i_trans];
      double ef = ef_arr[i_trans];
      if ((psi_i_prev == -1 && psi_f_prev == -1)) {
        strcpy(output_density_file, sp->out_file_base);
        sprintf(output_suffix, "_%d_%d.dens", psi_i, psi_f);
        strcat(output_density_file, output_suffix);
        out_file = fopen(output_density_file, "w+");
        fprintf(out_file, " Initial state #    %d E =  %8.5f 2xJ, 2xT =    %d   %d\n", psi_i + 1, ei, ji, ti);  
        fprintf(out_file, " Final state   #    %d E =  %8.5f 2xJ, 2xT =    %d   %d\n", psi_f + 1, ef, jf, tf);  

      } else if (psi_i != psi_i_prev || psi_f != psi_f_prev) {
        fprintf(out_file, "\n");
        fprintf(out_file, " Initial state #    %d E =  %8.5f 2xJ, 2xT =    %d   %d\n", psi_i + 1, ei, ji, ti);  
        fprintf(out_file, " Final state   #    %d E =  %8.5f 2xJ, 2xT =    %d   %d\n", psi_f + 1, ef, jf, tf);  

      }
      fprintf(out_file, " Jt =   %d, Tt = 0        1\n", j_op);
      int i_tcase;
      if (t_op == 0) {
        if (abs(ti + tf) != 0) {
	  i_tcase = 2;
          for (int i_orb1 = 0; i_orb1 < wd->n_orbits; i_orb1++) {
            for (int i_orb2 = 0; i_orb2 < wd->n_orbits; i_orb2++) {
              if (fabs(total[i_trans + sp->n_trans*(i_orb1 + i_orb2*wd->n_orbits)]) > pow(10, -12)) {
                fprintf(out_file, "    %d    %d  %8.5f  %8.5f\n", i_orb1 + 1, i_orb2 + 1, total[i_trans + sp->n_trans*(i_orb1 + i_orb2*wd->n_orbits)], total[i_trans + 1 + sp->n_trans*(i_orb1 + i_orb2*wd->n_orbits)]);
              }
            }
          }
	} else {
          i_tcase = 0;
          for (int i_orb1 = 0; i_orb1 < wd->n_orbits; i_orb1++) {
            for (int i_orb2 = 0; i_orb2 < wd->n_orbits; i_orb2++) {
              if (fabs(total[i_trans + sp->n_trans*(i_orb1 + i_orb2*wd->n_orbits)]) > pow(10, -12)) {
                fprintf(out_file, "    %d    %d  %8.5f  %8.5f\n", i_orb1 + 1, i_orb2 + 1, total[i_trans + sp->n_trans*(i_orb1 + i_orb2*wd->n_orbits)], 0.00000);
              }
            }
          }
        }
      } else {
	 i_tcase = 1;
         for (int i_orb1 = 0; i_orb1 < wd->n_orbits; i_orb1++) {
            for (int i_orb2 = 0; i_orb2 < wd->n_orbits; i_orb2++) {
              if (fabs(total[i_trans + sp->n_trans*(i_orb1 + i_orb2*wd->n_orbits)]) > pow(10, -12)) {
                fprintf(out_file, "    %d    %d  %8.5f  %8.5f\n", i_orb1 + 1, i_orb2 + 1, 0.00000, total[i_trans + sp->n_trans*(i_orb1 + i_orb2*wd->n_orbits)]);
              }
            }
          }
      }
      
      i_trans++;
      psi_i_prev = psi_i;
      psi_f_prev = psi_f;
      trans = trans->next;
      if (i_tcase == 2) {
        i_trans++;
	trans = trans->next;
      }
    }    
    fclose(out_file);
  }

  return;
}

void trace_two_body_nodes_dmtp1_13(int a, int b, int c, int d, int num_mj_1, float mj_min_1, int num_mj_2, float mj_min_2, int n_sds_int1, int* p1_array_f, sd_list** p2_list_i, sd_list** n1_list_f, wfnData* wd, int i_op, eigen_list* transition, double* density) {
  int ns = wd->n_shells;

  for (int imj1 = 0; imj1 < num_mj_1; imj1++) {
    float mj1 = imj1 + mj_min_1;
    for (int imj2 = 0; imj2 < num_mj_2; imj2++) {
      float mj2 = imj2 + mj_min_2;
      if ((mj1 + mj2 != 0) && (mj1 + mj2 != 0.5)) {continue;}
      for (int ipar1 = 0; ipar1 <= 1; ipar1++) {
        int ipar2;
	if (wd->parity_i == '+') {
	  if (ipar1 == 0) {
		  ipar2 = 0;
	  } else {
		  ipar2 = 1;
	  }
	} else if (wd->parity_i == '-') {
		if (ipar1 == 0) {
			ipar2 = 1;
		} else {
			ipar2 = 0;
		}
	} else {printf("Parity error\n"); exit(0);}

        sd_list* node1 = p2_list_i[ipar1 + 2*(imj1 + num_mj_1*(d + c*ns))];
        while (node1 != NULL) {
          unsigned int ppn = node1->pn;
          unsigned int ppi = node1->pi;
          int phase1 = node1->phase;
          int ppf = p1_array_f[(ppn - 1) + n_sds_int1*b];
          if (ppf == 0) {node1 = node1->next; continue;}
          int phase2 = 1;
          node1 = node1->next;
          if (ppf < 0) {
            ppf *= -1;
            phase2 = -1;
          }
          sd_list* node2 = n1_list_f[ipar2 + 2*(imj2 + num_mj_2*a)];
          while (node2 != NULL) {
            int pnf = node2->pn;
	    int pni = node2->pi;
	    int phase3 = node2->phase;
            int index_i = -1;
            int index_f = -1;
            if (i_op == 0) {
              unsigned int p_hash_i = ppi + wd->n_sds_p_i*(pni % HASH_SIZE);
              wh_list* node3 = wd->wh_hash_i[p_hash_i];
              while (node3 != NULL) {
                if ((pni == node3->pn) && (ppi == node3->pp)) {
                  index_i = node3->index;
                  break;
                }
                node3 = node3->next;
              }
              if (index_i < 0) {node2 = node2->next; continue;}
              unsigned int p_hash_f = ppf + wd->n_sds_p_f*(pnf % HASH_SIZE);
              node3 = wd->wh_hash_f[p_hash_f];
              while (node3 != NULL) {
                if ((pnf == node3->pn) && (ppf == node3->pp)) {
                  index_f = node3->index;
                  break;
                }
                node3 = node3->next;
              }
              if (index_f < 0) {node2 = node2->next; continue;}
              eigen_list* eig_pair = transition;
              int i_trans = 0;
              while (eig_pair != NULL) {
                int psi_i = eig_pair->eig_i;
                int psi_f = eig_pair->eig_f;
	        density[i_trans] += wd->bc_i[psi_i + wd->n_eig_i*index_i]*wd->bc_f[psi_f + wd->n_eig_f*index_f]*phase1*phase2*phase3;
                i_trans++;
                eig_pair = eig_pair->next;
              }
             } else {
              unsigned int p_hash_i = pni + wd->n_sds_p_i*(ppi % HASH_SIZE);
              wh_list* node3 = wd->wh_hash_i[p_hash_i];
              while (node3 != NULL) {
                if ((pni == node3->pp) && (ppi == node3->pn)) {
                  index_i = node3->index;
                  break;
                }
                node3 = node3->next;
              }
              if (index_i < 0) {node2 = node2->next; continue;}

              unsigned int p_hash_f = pnf + wd->n_sds_p_f*(ppf % HASH_SIZE);
              node3 = wd->wh_hash_f[p_hash_f];
              while (node3 != NULL) {
                if ((pnf == node3->pp) && (ppf == node3->pn)) {
                  index_f = node3->index;
                  break;
                }
                node3 = node3->next;
              }
              if (index_f < 0) {node2 = node2->next; continue;}
              eigen_list* eig_pair = transition;
              int i_trans = 0;
              while (eig_pair != NULL) {
                int psi_i = eig_pair->eig_i;
                int psi_f = eig_pair->eig_f;
	        density[i_trans] += wd->bc_i[psi_i + wd->n_eig_i*index_i]*wd->bc_f[psi_f + wd->n_eig_f*index_f]*phase1*phase2*phase3;
                i_trans++;
                eig_pair = eig_pair->next;
              }
            } 
            node2 = node2->next;
	  }
	} 
      }
    }
   }
  return;
}


void trace_two_body_nodes_dmtp1_31(int a, int b, int c, int d, int num_mj_1, float mj_min_1, int num_mj_2, float mj_min_2, int n_sds_int2, int* p2_array_f, sd_list** p1_list_i, sd_list** n1_list_i, wfnData* wd, int i_op, eigen_list* transition, double* density) {
  int ns = wd->n_shells;
 
  for (int imj1 = 0; imj1 < num_mj_1; imj1++) {
    float mj1 = imj1 + mj_min_1;
    for (int imj2 = 0; imj2 < num_mj_2; imj2++) {
      float mj2 = imj2 + mj_min_2;
      if ((mj1 + mj2 != 0) && (mj1 + mj2 != 0.5)) {continue;}
      for (int ipar1 = 0; ipar1 <= 1; ipar1++) {
        int ipar2;
	if (wd->parity_i == '+') {
	  if (ipar1 == 0) {
		  ipar2 = 0;
	  } else {
		  ipar2 = 1;
	  }
	} else if (wd->parity_i == '-') {
		if (ipar1 == 0) {
			ipar2 = 1;
		} else {
			ipar2 = 0;
		}
	} else {printf("Parity error\n"); exit(0);}
        
	sd_list* node1 = p1_list_i[ipar1 + 2*(imj1 + num_mj_1*c)];
        while (node1 != NULL) {
          unsigned int ppn = node1->pn;
          unsigned int ppi = node1->pi;
          int phase1 = node1->phase;
          int ppf = p2_array_f[(ppn - 1) + n_sds_int2*(a + ns*b)];
          if (ppf == 0) {node1 = node1->next; continue;}
          int phase2 = 1;
          node1 = node1->next;
          if (ppf < 0) {
            ppf *= -1;
            phase2 = -1;
          }
          sd_list* node2 = n1_list_i[ipar2 + 2*(imj2 + num_mj_2*d)];

          while (node2 != NULL) {
            int pnf = node2->pn;
	    int pni = node2->pi;
	    int phase3 = node2->phase;
            int index_i = -1;
            int index_f = -1;
            if (i_op == 0) {
              unsigned int p_hash_i = ppi + wd->n_sds_p_i*(pni % HASH_SIZE);
              wh_list* node3 = wd->wh_hash_i[p_hash_i];
              while (node3 != NULL) {
                if ((pni == node3->pn) && (ppi == node3->pp)) {
                  index_i = node3->index;
                  break;
                }
                node3 = node3->next;
              }
              if (index_i < 0) {node2 = node2->next; continue;}
              unsigned int p_hash_f = ppf + wd->n_sds_p_f*(pnf % HASH_SIZE);
              node3 = wd->wh_hash_f[p_hash_f];
              while (node3 != NULL) {
                if ((pnf == node3->pn) && (ppf == node3->pp)) {
                  index_f = node3->index;
                  break;
                }
                node3 = node3->next;
              }
              if (index_f < 0) {node2 = node2->next; continue;}
              eigen_list* eig_pair = transition;
              int i_trans = 0;
              while (eig_pair != NULL) {
                int psi_i = eig_pair->eig_i;
                int psi_f = eig_pair->eig_f;
	        density[i_trans] += wd->bc_i[psi_i + wd->n_eig_i*index_i]*wd->bc_f[psi_f + wd->n_eig_f*index_f]*phase1*phase2*phase3;
                i_trans++;
                eig_pair = eig_pair->next;
              }
             } else {
              unsigned int p_hash_i = pni + wd->n_sds_p_i*(ppi % HASH_SIZE);
              wh_list* node3 = wd->wh_hash_i[p_hash_i];
              while (node3 != NULL) {
                if ((pni == node3->pp) && (ppi == node3->pn)) {
                  index_i = node3->index;
                  break;
                }
                node3 = node3->next;
              }
              if (index_i < 0) {node2 = node2->next; continue;}

              unsigned int p_hash_f = pnf + wd->n_sds_p_f*(ppf % HASH_SIZE);
              node3 = wd->wh_hash_f[p_hash_f];
              while (node3 != NULL) {
                if ((pnf == node3->pp) && (ppf == node3->pn)) {
                  index_f = node3->index;
                  break;
                }
                node3 = node3->next;
              }
              if (index_f < 0) {node2 = node2->next; continue;}
              eigen_list* eig_pair = transition;
              int i_trans = 0;
              while (eig_pair != NULL) {
                int psi_i = eig_pair->eig_i;
                int psi_f = eig_pair->eig_f;
	        density[i_trans] += wd->bc_i[psi_i + wd->n_eig_i*index_i]*wd->bc_f[psi_f + wd->n_eig_f*index_f]*phase1*phase2*phase3;
                i_trans++;
                eig_pair = eig_pair->next;
              }
            } 
            node2 = node2->next;
        
	  }
	} 
      }
    }
   }
  return;
}


void trace_two_body_nodes_dmt0_40(int a, int b, int c, int d, int num_mj_p_i, float mj_min_p_i, int num_mj_n_i, float mj_min_n_i, int n_sds_int2, int* p2_array_f, sd_list** p2_list_i, wf_list** n0_list_i, wfnData* wd, int i_op, eigen_list* transition, double* density) {
  int ns = wd->n_shells;


  for (int imjp = 0; imjp < num_mj_p_i; imjp++) {
    float mjp = imjp + mj_min_p_i;
    for (int imjn = 0; imjn < num_mj_n_i; imjn++) {
      float mjn = imjn + mj_min_n_i;
      if ((mjp + mjn != 0) && (mjp + mjn != 0.5)) {continue;}
      for (int ipar_p = 0; ipar_p <= 1; ipar_p++) {
        int ipar_n;
	if (wd->parity_i == '+') {
	  if (ipar_p == 0) {
		  ipar_n = 0;
	  } else {
		  ipar_n = 1;
	  }
	} else if (wd->parity_i == '-') {
		if (ipar_p == 0) {
			ipar_n = 1;
		} else {
			ipar_n = 0;
		}
	} else {printf("Parity error\n"); exit(0);}


        sd_list* node1 = p2_list_i[ipar_p + 2*(imjp + num_mj_p_i*(c + d*ns))];
        while (node1 != NULL) {
          unsigned int ppn = node1->pn;
          unsigned int ppi = node1->pi;
          int phase1 = node1->phase;
          int ppf = p2_array_f[(ppn - 1) + n_sds_int2*(a + ns*b)];
          if (ppf == 0) {node1 = node1->next; continue;}
          int phase2 = 1;
          node1 = node1->next;
          if (ppf < 0) {
            ppf *= -1;
            phase2 = -1;
          }
          wf_list* node2 = n0_list_i[ipar_n + 2*imjn];
          while (node2 != NULL) {
            int pn = node2->p;
            int index_i = -1;
            int index_f = -1;

            if (i_op == 0) {
              unsigned int p_hash_i = ppi + wd->n_sds_p_i*(pn % HASH_SIZE);
              wh_list* node3 = wd->wh_hash_i[p_hash_i];
              while (node3 != NULL) {
                if ((pn == node3->pn) && (ppi == node3->pp)) {
                  index_i = node3->index;
                  break;
                }
                node3 = node3->next;
              }
              if (index_i < 0) {node2 = node2->next; continue;}

              unsigned int p_hash_f = ppf + wd->n_sds_p_f*(pn % HASH_SIZE);
              node3 = wd->wh_hash_f[p_hash_f];
              while (node3 != NULL) {
                if ((pn == node3->pn) && (ppf == node3->pp)) {
                  index_f = node3->index;
                  break;
                }
                node3 = node3->next;
              }
              if (index_f < 0) {node2 = node2->next; continue;}
              eigen_list* eig_pair = transition;
              int i_trans = 0;
              while (eig_pair != NULL) {
                int psi_i = eig_pair->eig_i;
                int psi_f = eig_pair->eig_f;
	        density[i_trans] += wd->bc_i[psi_i + wd->n_eig_i*index_i]*wd->bc_f[psi_f + wd->n_eig_f*index_f]*phase1*phase2;
                i_trans++;
                eig_pair = eig_pair->next;
              }
             } else {
              unsigned int p_hash_i = pn + wd->n_sds_p_i*(ppi % HASH_SIZE);
              wh_list* node3 = wd->wh_hash_i[p_hash_i];
              while (node3 != NULL) {
                if ((pn == node3->pp) && (ppi == node3->pn)) {
                  index_i = node3->index;
                  break;
                }
                node3 = node3->next;
              }
              if (index_i < 0) {node2 = node2->next; continue;}

              unsigned int p_hash_f = pn + wd->n_sds_p_f*(ppf % HASH_SIZE);
              node3 = wd->wh_hash_f[p_hash_f];
              while (node3 != NULL) {
                if ((pn == node3->pp) && (ppf == node3->pn)) {
                  index_f = node3->index;
                  break;
                }
                node3 = node3->next;
              }
              if (index_f < 0) {node2 = node2->next; continue;}
              eigen_list* eig_pair = transition;
              int i_trans = 0;
              while (eig_pair != NULL) {
                int psi_i = eig_pair->eig_i;
                int psi_f = eig_pair->eig_f;
	        density[i_trans] += wd->bc_i[psi_i + wd->n_eig_i*index_i]*wd->bc_f[psi_f + wd->n_eig_f*index_f]*phase1*phase2;
                i_trans++;
                eig_pair = eig_pair->next;
              }
            } 
            node2 = node2->next;
          }
	} 
      }
    }
   }
  return;
}

void trace_two_body_nodes_dmtp2(int a, int b, int c, int d, int num_mj_p_i, float mj_min_p_i, int num_mj_n_i, float mj_min_n_i, sd_list** n2_list_i, sd_list** p2_list_f, wfnData* wd, int i_op, eigen_list* transition, double* density) {
/* 

*/
  int ns = wd->n_shells;
  for (int imjp = 0; imjp < num_mj_p_i; imjp++) {
    float mjp = imjp + mj_min_p_i;
    for (int imjn = 0; imjn < num_mj_n_i; imjn++) {
      float mjn = imjn + mj_min_n_i;
      if ((mjp + mjn != 0) && (mjp + mjn != 0.5)) {continue;}
      for (int ipar_p = 0; ipar_p <= 1; ipar_p++) {
        int ipar_n;
	if (wd->parity_i == '+') {
	  if (ipar_p == 0) {
		  ipar_n = 0;
	  } else {
		  ipar_n = 1;
	  }
	} else if (wd->parity_i == '-') {
		if (ipar_p == 0) {
			ipar_n = 1;
		} else {
			ipar_n = 0;
		}
	} else {printf("Parity error\n"); exit(0);}

        sd_list* node1 = n2_list_i[ipar_n + 2*(imjn + num_mj_n_i*(c + d*ns))];
        // Loop over final states resulting from 2x a_op
        while (node1 != NULL) {
          unsigned int pnf = node1->pn; // Get final state p_f
          unsigned int pni = node1->pi; // Get initial state p_i
          int phase1 = node1->phase;
          // Get list of p_i associated to p_f
          sd_list* node2 = p2_list_f[ipar_p + 2*(imjp + num_mj_p_i*(a + b*ns))]; //hash corresponds to a_op operators
          // Loop over n_f  
          while (node2 != NULL) {
            unsigned int ppf = node2->pi;
            unsigned int ppi = node2->pn;
            int phase2 = node2->phase;
            int index_i = -1;
            int index_f = -1;
            if (i_op == 0) {
              unsigned int p_hash_i = ppi + wd->n_sds_p_i*(pni % HASH_SIZE);
              wh_list* node3 = wd->wh_hash_i[p_hash_i];
              while (node3 != NULL) {
                if ((pni == node3->pn) && (ppi == node3->pp)) {
                  index_i = node3->index;
                  break;
                }
                node3 = node3->next;
              }
              if (index_i < 0) {node2 = node2->next; continue;}

              unsigned int p_hash_f = ppf + wd->n_sds_p_f*(pnf % HASH_SIZE);
              node3 = wd->wh_hash_f[p_hash_f];
              while (node3 != NULL) {
                if ((pnf == node3->pn) && (ppf == node3->pp)) {
                  index_f = node3->index;
                  break;
                }
                node3 = node3->next;
              }
              if (index_f < 0) {node2 = node2->next; continue;}
              eigen_list* eig_pair = transition;
              int i_trans = 0;
              while (eig_pair != NULL) {
                int psi_i = eig_pair->eig_i;
                int psi_f = eig_pair->eig_f;
	        density[i_trans] += wd->bc_i[psi_i + wd->n_eig_i*index_i]*wd->bc_f[psi_f + wd->n_eig_f*index_f]*phase1*phase2;
                i_trans++;
                eig_pair = eig_pair->next;
              }
             } else {
              unsigned int p_hash_i = pni + wd->n_sds_p_i*(ppi % HASH_SIZE);
              wh_list* node3 = wd->wh_hash_i[p_hash_i];
              while (node3 != NULL) {
                if ((pni == node3->pp) && (ppi == node3->pn)) {
                  index_i = node3->index;
                  break;
                }
                node3 = node3->next;
              }
              if (index_i < 0) {node2 = node2->next; continue;}

              unsigned int p_hash_f = pnf + wd->n_sds_p_f*(ppf % HASH_SIZE);
              node3 = wd->wh_hash_f[p_hash_f];
              while (node3 != NULL) {
                if ((pnf == node3->pp) && (ppf == node3->pn)) {
                  index_f = node3->index;
                  break;
                }
                node3 = node3->next;
              }
              if (index_f < 0) {node2 = node2->next; continue;}
              eigen_list* eig_pair = transition;
              int i_trans = 0;
              while (eig_pair != NULL) {
                int psi_i = eig_pair->eig_i;
                int psi_f = eig_pair->eig_f;
	        density[i_trans] += wd->bc_i[psi_i + wd->n_eig_i*index_i]*wd->bc_f[psi_f + wd->n_eig_f*index_f]*phase1*phase2;
                i_trans++;
                eig_pair = eig_pair->next;
              }
            } 
            node2 = node2->next;
          }
          node1 = node1->next;
        }
      }
    }  
  }
  return;
}

void trace_two_body_nodes_dmt0_22_alt(int a, int b, int c, int d, int num_mj_p_i, float mj_min_p_i, float mj_min_p_f, float mj_max_p_f, int num_mj_n_i, float mj_min_n_i, float mj_min_n_f, float mj_max_n_f, sd_list** p2_list_if, sd_list** n2_list_if, wfnData* wd, eigen_list* transition, double* density, int* jz_shell) {
 
  int ns = wd->n_shells;
  for (int imjp = 0; imjp < num_mj_p_i; imjp++) {
    float mjp = imjp + mj_min_p_i;
    if (mjp + jz_shell[a]/2.0 - jz_shell[b]/2.0 > mj_max_p_f || mjp + jz_shell[a]/2.0 - jz_shell[b]/2.0 < mj_min_p_f) {continue;}
    for (int imjn = 0; imjn < num_mj_n_i; imjn++) {
      float mjn = imjn + mj_min_n_i;
      if ((mjp + mjn != 0) && (mjp + mjn != 0.5)) {continue;}
      if (mjn + jz_shell[c]/2.0 - jz_shell[d]/2.0 > mj_max_n_f || mjn + jz_shell[c]/2.0 - jz_shell[d]/2.0 < mj_min_n_f) {continue;}
      for (int ipar_p = 0; ipar_p <= 1; ipar_p++) {
        int ipar_n;
	if (wd->parity_i == '+') {
	  if (ipar_p == 0) {
		  ipar_n = 0;
	  } else {
		  ipar_n = 1;
	  }
	} else if (wd->parity_i == '-') {
		if (ipar_p == 0) {
			ipar_n = 1;
		} else {
			ipar_n = 0;
		}
	} else {printf("Parity error\n"); exit(0);}

        sd_list* node_pi = p2_list_if[ipar_p + 2*(imjp + num_mj_p_i*(b + ns*a))]; // Get proton states resulting from p_b |p_i>
        while (node_pi != NULL) {
          int ppi = node_pi->pi; 
          int ppf = node_pi->pn; // Get pn = p_a |p_i>
          int phase1 = node_pi->phase;
          
	  sd_list* node_ni = n2_list_if[ipar_n + 2*(imjn + num_mj_n_i*(d + ns*c))]; // Get neutron states resulting from n_d |n_i>
          while (node_ni != NULL) {
            int pni = node_ni->pi;
            int pnf = node_ni->pn; // Get nn = n_d |n_i>
            int phase2 = node_ni->phase;
            
	    unsigned int p_hash_i = ppi + wd->n_sds_p_i*(pni % HASH_SIZE);
            unsigned int p_hash_f = ppf + wd->n_sds_p_f*(pnf % HASH_SIZE);
            wh_list *node3 = wd->wh_hash_i[p_hash_i]; // hash corresponds to a_op operators
   
            int index_i = -1;
            int index_f = -1;
            while (node3 != NULL) {
              if ((pni == node3->pn) && (ppi == node3->pp)) {
                index_i = node3->index;
                break;
              }
              node3 = node3->next;
            }
            if (index_i < 0) {node_ni = node_ni->next; continue;}
            node3 = wd->wh_hash_f[p_hash_f];
            while (node3 != NULL) {
              if ((pnf == node3->pn) && (ppf == node3->pp)) {
                index_f = node3->index;
                break;
              }
              node3 = node3->next;
            }
            if (index_f < 0) {node_ni = node_ni->next; continue;}
            eigen_list* eig_pair = transition;
            int i_trans = 0;
            while (eig_pair != NULL) {
              int psi_i = eig_pair->eig_i;
              int psi_f = eig_pair->eig_f;
              density[i_trans] += wd->bc_i[psi_i + wd->n_eig_i*index_i]*wd->bc_f[psi_f + wd->n_eig_f*index_f]*phase1*phase2;
              i_trans++;
              eig_pair = eig_pair->next;
            }
            node_ni = node_ni->next;
          }
          node_pi = node_pi->next;
        }
      }
    } 
  } 
  return;
}


void trace_two_body_nodes_dmt0_22(int a, int b, int c, int d, int num_mj_p_i, float mj_min_p_i, float mj_min_p_f, float mj_max_p_f, int num_mj_n_i, float mj_min_n_i, float mj_min_n_f, float mj_max_n_f, int n_sds_p_int1, int n_sds_n_int1, sd_list** p1_list_i, sd_list** n1_list_i, int* p1_array_f, int* n1_array_f, wfnData* wd, eigen_list* transition, double* density, int* jz_shell) {
  
  for (int imjp = 0; imjp < num_mj_p_i; imjp++) {
    float mjp = imjp + mj_min_p_i;
    if (mjp + jz_shell[a]/2.0 - jz_shell[b]/2.0 > mj_max_p_f || mjp + jz_shell[a]/2.0 - jz_shell[b]/2.0 < mj_min_p_f) {continue;}
    for (int imjn = 0; imjn < num_mj_n_i; imjn++) {
      float mjn = imjn + mj_min_n_i;
      if ((mjp + mjn != 0) && (mjp + mjn != 0.5)) {continue;}
      if (mjn + jz_shell[c]/2.0 - jz_shell[d]/2.0 > mj_max_n_f || mjn + jz_shell[c]/2.0 - jz_shell[d]/2.0 < mj_min_n_f) {continue;}
      for (int ipar_p = 0; ipar_p <= 1; ipar_p++) {
        int ipar_n;
	if (wd->parity_i == '+') {
	  if (ipar_p == 0) {
		  ipar_n = 0;
	  } else {
		  ipar_n = 1;
	  }
	} else if (wd->parity_i == '-') {
		if (ipar_p == 0) {
			ipar_n = 1;
		} else {
			ipar_n = 0;
		}
	} else {printf("Parity error\n"); exit(0);}

        sd_list* node_pi = p1_list_i[ipar_p + 2*(imjp + num_mj_p_i*b)]; // Get proton states resulting from p_b |p_i>
        while (node_pi != NULL) {
          int ppi = node_pi->pi; 
          int ppn = node_pi->pn; // Get pn = p_a |p_i>
          int phase1 = node_pi->phase;
          int phase2 = 1;
          int ppf = p1_array_f[(ppn - 1) + n_sds_p_int1*a]; // Get pf such that pn = p_a |p_f>
          if (ppf == 0) {node_pi = node_pi->next; continue;}
          if (ppf < 0) {
            ppf *= -1;
            phase2 = -1;
          }
          sd_list* node_ni = n1_list_i[ipar_n + 2*(imjn + num_mj_n_i*d)]; // Get neutron states resulting from n_d |n_i>
          while (node_ni != NULL) {
            int pni = node_ni->pi;
            int pnn = node_ni->pn; // Get nn = n_d |n_i>
            int phase3 = node_ni->phase;
            int phase4 = 1;
            int pnf = n1_array_f[(pnn - 1) + n_sds_n_int1*c];
            if (pnf == 0) {node_ni = node_ni->next; continue;}
            if (pnf < 0) {
              pnf *= -1;
              phase4 = -1;
            }
            unsigned int p_hash_i = ppi + wd->n_sds_p_i*(pni % HASH_SIZE);
            unsigned int p_hash_f = ppf + wd->n_sds_p_f*(pnf % HASH_SIZE);
            wh_list *node3 = wd->wh_hash_i[p_hash_i]; // hash corresponds to a_op operators
   
            int index_i = -1;
            int index_f = -1;
            while (node3 != NULL) {
              if ((pni == node3->pn) && (ppi == node3->pp)) {
                index_i = node3->index;
                break;
              }
              node3 = node3->next;
            }
            if (index_i < 0) {node_ni = node_ni->next; continue;}
            node3 = wd->wh_hash_f[p_hash_f];
            while (node3 != NULL) {
              if ((pnf == node3->pn) && (ppf == node3->pp)) {
                index_f = node3->index;
                break;
              }
              node3 = node3->next;
            }
            if (index_f < 0) {node_ni = node_ni->next; continue;}
            eigen_list* eig_pair = transition;
            int i_trans = 0;
            while (eig_pair != NULL) {
              int psi_i = eig_pair->eig_i;
              int psi_f = eig_pair->eig_f;
              density[i_trans] += wd->bc_i[psi_i + wd->n_eig_i*index_i]*wd->bc_f[psi_f + wd->n_eig_f*index_f]*phase1*phase2*phase3*phase4;
              i_trans++;
              eig_pair = eig_pair->next;
            }
            node_ni = node_ni->next;
          }
          node_pi = node_pi->next;
        }
      }
    } 
  } 
  return;
}

void build_two_body_jumps_dmt0_alt(int n_s, int n_p, float mj_min, float mj_max, int num_mj, int n_sds_i, int n_sds_int1, int n_sds_int2, int*a1_array_f, int* a2_array_f, wf_list** a0_list_i, sd_list** a1_list_i, sd_list** a2_list_i, sd_list** a2_list_if, int* jz_shell, int* l_shell) {
/*
  Input(s):
    mj_min_i: 
*/
  for (int j = 1; j <= n_sds_i; j++) {
    float mj = m_from_p(j, n_s, n_p, jz_shell);
    if ((mj < mj_min) || (mj > mj_max)) {continue;}
    int i_parity = (parity_from_p(j, n_s, n_p, l_shell) + 1)/2;
    int i_mj = mj - mj_min;
    if (a0_list_i[i_parity + 2*i_mj] == NULL) {
      a0_list_i[i_parity + 2*i_mj] = create_wf_node(j, NULL);
    } else {
      wf_append(a0_list_i[i_parity + 2*i_mj], j);
    }
    int j_min = j_min_from_p(n_s, n_p, j);
    for (int b = j_min - 1; b < n_s; b++) {
      int phase1;
      int pn1 = a_op(n_s, n_p, j, b + 1, &phase1, j_min);
      if (pn1 == 0) {continue;}
      for (int a = j_min - 1; a < b; a++) {
        int phase2;
        int pn2 = a_op(n_s, n_p - 1, pn1, a + 1, &phase2, j_min);
        if (pn2 == 0) {continue;}
        if (a2_list_i[i_parity + 2*(i_mj + num_mj*(b + a*n_s))] == NULL) {
          a2_list_i[i_parity + 2*(i_mj + num_mj*(b + a*n_s))] = create_sd_node(j, pn2, phase1*phase2, NULL);
        } else {
          sd_append(a2_list_i[i_parity + 2*(i_mj + num_mj*(b + a*n_s))], j, pn2, phase1*phase2);
        }
        if (a2_list_i[i_parity + 2*(i_mj + num_mj*(a + b*n_s))] == NULL) {
          a2_list_i[i_parity + 2*(i_mj + num_mj*(a + b*n_s))] = create_sd_node(j, pn2, -phase1*phase2, NULL);
        } else {
          sd_append(a2_list_i[i_parity + 2*(i_mj + num_mj*(a + b*n_s))], j, pn2, -phase1*phase2);
        }
        a2_array_f[(pn2 - 1) + n_sds_int2*(b + a*n_s)] = phase1*phase2*j;
        a2_array_f[(pn2 - 1) + n_sds_int2*(a + b*n_s)] = -phase1*phase2*j;
      }
      for (int a = 0; a < n_s; a++) {
        int phase2;
        int pn2 = a_op_dag(n_s, n_p - 1, pn1, a + 1, &phase2, 1);
        if (pn2 == 0) {continue;}
        float mj = m_from_p(pn2, n_s, n_p, jz_shell);
        if ((mj < mj_min) || (mj > mj_max)) {continue;}

        if (a2_list_if[i_parity + 2*(i_mj + num_mj*(b + a*n_s))] == NULL) {
          a2_list_if[i_parity + 2*(i_mj + num_mj*(b + a*n_s))] = create_sd_node(j, pn2, phase1*phase2, NULL);
        } else {
          sd_append(a2_list_if[i_parity + 2*(i_mj + num_mj*(b + a*n_s))], j, pn2, phase1*phase2);
        }
      }
    }
  } 

  return;
}


void build_two_body_jumps_dmt0(int n_s, int n_p, float mj_min, float mj_max, int num_mj, int n_sds_i, int n_sds_int1, int n_sds_int2, int*a1_array_f, int* a2_array_f, wf_list** a0_list_i, sd_list** a1_list_i, sd_list** a2_list_i, int* jz_shell, int* l_shell) {
/*
  Input(s):
    mj_min_i: 
*/
  for (int j = 1; j <= n_sds_i; j++) {
    float mj = m_from_p(j, n_s, n_p, jz_shell);
    if ((mj < mj_min) || (mj > mj_max)) {continue;}
    int i_parity = (parity_from_p(j, n_s, n_p, l_shell) + 1)/2;
    int i_mj = mj - mj_min;
    if (a0_list_i[i_parity + 2*i_mj] == NULL) {
      a0_list_i[i_parity + 2*i_mj] = create_wf_node(j, NULL);
    } else {
      wf_append(a0_list_i[i_parity + 2*i_mj], j);
    }
    int j_min = j_min_from_p(n_s, n_p, j);
    for (int b = j_min - 1; b < n_s; b++) {
      int phase1;
      int pn1 = a_op(n_s, n_p, j, b + 1, &phase1, j_min);
      if (pn1 == 0) {continue;}
      a1_array_f[(pn1 - 1) + b*n_sds_int1] = j*phase1;
      if (a1_list_i[i_parity + 2*(i_mj + num_mj*b)] == NULL) {
        a1_list_i[i_parity + 2*(i_mj + num_mj*b)] = create_sd_node(j, pn1, phase1, NULL);
      } else {
        sd_append(a1_list_i[i_parity + 2*(i_mj + num_mj*b)], j, pn1, phase1);
      }
      for (int a = j_min - 1; a < b; a++) {
        int phase2;
        int pn2 = a_op(n_s, n_p - 1, pn1, a + 1, &phase2, j_min);
        if (pn2 == 0) {continue;}
        if (a2_list_i[i_parity + 2*(i_mj + num_mj*(b + a*n_s))] == NULL) {
          a2_list_i[i_parity + 2*(i_mj + num_mj*(b + a*n_s))] = create_sd_node(j, pn2, phase1*phase2, NULL);
        } else {
          sd_append(a2_list_i[i_parity + 2*(i_mj + num_mj*(b + a*n_s))], j, pn2, phase1*phase2);
        }
        if (a2_list_i[i_parity + 2*(i_mj + num_mj*(a + b*n_s))] == NULL) {
          a2_list_i[i_parity + 2*(i_mj + num_mj*(a + b*n_s))] = create_sd_node(j, pn2, -phase1*phase2, NULL);
        } else {
          sd_append(a2_list_i[i_parity + 2*(i_mj + num_mj*(a + b*n_s))], j, pn2, -phase1*phase2);
        }
        a2_array_f[(pn2 - 1) + n_sds_int2*(b + a*n_s)] = phase1*phase2*j;
        a2_array_f[(pn2 - 1) + n_sds_int2*(a + b*n_s)] = -phase1*phase2*j;
      }
    }
  } 

  return;
}


void build_two_body_jumps_dmt0_trunc(int n_s, int n_p, float mj_min, float mj_max, int num_mj, int n_sds_i, int n_sds_int1, int n_sds_int2, int*a1_array_f, int* a2_array_f, wf_list** a0_list_i, sd_list** a1_list_i, sd_list** a2_list_i, int* jz_shell, int* l_shell, int* w_shell, int w_max) {
/*
  Input(s):
    mj_min_i: 
*/
  for (int j = 1; j <= n_sds_i; j++) {
    if (w_from_p(j, n_s, n_p, w_shell) > w_max) {continue;}
    float mj = m_from_p(j, n_s, n_p, jz_shell);
    if ((mj < mj_min) || (mj > mj_max)) {continue;}
    int i_parity = (parity_from_p(j, n_s, n_p, l_shell) + 1)/2;
    int i_mj = mj - mj_min;
    if (a0_list_i[i_parity + 2*i_mj] == NULL) {
      a0_list_i[i_parity + 2*i_mj] = create_wf_node(j, NULL);
    } else {
      wf_append(a0_list_i[i_parity + 2*i_mj], j);
    }
    int j_min = j_min_from_p(n_s, n_p, j);
    for (int b = j_min - 1; b < n_s; b++) {
      int phase1;
      int pn1 = a_op(n_s, n_p, j, b + 1, &phase1, j_min);
      if (pn1 == 0) {continue;}
      a1_array_f[(pn1 - 1) + b*n_sds_int1] = j*phase1;
      if (a1_list_i[i_parity + 2*(i_mj + num_mj*b)] == NULL) {
        a1_list_i[i_parity + 2*(i_mj + num_mj*b)] = create_sd_node(j, pn1, phase1, NULL);
      } else {
        sd_append(a1_list_i[i_parity + 2*(i_mj + num_mj*b)], j, pn1, phase1);
      }
      for (int a = j_min - 1; a < b; a++) {
        int phase2;
        int pn2 = a_op(n_s, n_p - 1, pn1, a + 1, &phase2, j_min);
        if (pn2 == 0) {continue;}
        if (a2_list_i[i_parity + 2*(i_mj + num_mj*(b + a*n_s))] == NULL) {
          a2_list_i[i_parity + 2*(i_mj + num_mj*(b + a*n_s))] = create_sd_node(j, pn2, phase1*phase2, NULL);
        } else {
          sd_append(a2_list_i[i_parity + 2*(i_mj + num_mj*(b + a*n_s))], j, pn2, phase1*phase2);
        }
        if (a2_list_i[i_parity + 2*(i_mj + num_mj*(a + b*n_s))] == NULL) {
          a2_list_i[i_parity + 2*(i_mj + num_mj*(a + b*n_s))] = create_sd_node(j, pn2, -phase1*phase2, NULL);
        } else {
          sd_append(a2_list_i[i_parity + 2*(i_mj + num_mj*(a + b*n_s))], j, pn2, -phase1*phase2);
        }
        a2_array_f[(pn2 - 1) + n_sds_int2*(b + a*n_s)] = phase1*phase2*j;
        a2_array_f[(pn2 - 1) + n_sds_int2*(a + b*n_s)] = -phase1*phase2*j;
      }
    }
  } 

  return;
}


void build_two_body_jumps_dmtp2_trunc(int n_s, int n_proton, int num_mj_p_i, float mj_min_p_i, float mj_max_p_i, int num_mj_p_f, float mj_min_p_f, float mj_max_p_f, int n_sds_p_f, sd_list** p2_list_f, int n_neutron, int num_mj_n_i, float mj_min_n_i, float mj_max_n_i, int num_mj_n_f, float mj_min_n_f, float mj_max_n_f, int n_sds_n_i, sd_list** n2_list_i, int* jz_shell, int* l_shell, int* w_shell, int w_max_p_i, int w_max_p_f, int w_max_n_i, int w_max_n_f) {

  for (int j = 1; j <= n_sds_n_i; j++) {
    if (w_from_p(j, n_s, n_neutron, w_shell) > w_max_n_i) {continue;}
    float mj = m_from_p(j, n_s, n_neutron, jz_shell);
    if ((mj < mj_min_n_i) || (mj > mj_max_n_i)) {continue;}
    int i_mj = mj - mj_min_n_i;
    int i_parity = (parity_from_p(j, n_s, n_neutron, l_shell) + 1)/2;
    int j_min = j_min_from_p(n_s, n_neutron, j);
    
    for (int b = j_min - 1; b < n_s; b++) {
      int phase1;
      int pn1 = a_op(n_s, n_neutron, j, b + 1, &phase1, j_min);
      if (pn1 == 0) {continue;}
      for (int a = j_min - 1; a < b; a++) {
        int phase2;
        int pn2 = a_op(n_s, n_neutron - 1, pn1, a + 1, &phase2, j_min);
        if (pn2 == 0) {continue;}
        mj = m_from_p(pn2, n_s, n_neutron - 2, jz_shell);
        if ((mj < mj_min_n_f) || (mj > mj_max_n_f)) {continue;}
        if (w_from_p(pn2, n_s, n_neutron - 2, w_shell) > w_max_n_f) {continue;}

        if (n2_list_i[i_parity + 2*(i_mj + num_mj_n_i*(b + a*n_s))] == NULL) {
          n2_list_i[i_parity + 2*(i_mj + num_mj_n_i*(b + a*n_s))] = create_sd_node(j, pn2, phase1*phase2, NULL);
        } else {
          sd_append(n2_list_i[i_parity + 2*(i_mj + num_mj_n_i*(b + a*n_s))], j, pn2, phase1*phase2);
        }
        if (n2_list_i[i_parity + 2*(i_mj + num_mj_n_i*(a + b*n_s))] == NULL) {
          n2_list_i[i_parity + 2*(i_mj + num_mj_n_i*(a + b*n_s))] = create_sd_node(j, pn2, -phase1*phase2, NULL);
        } else {
          sd_append(n2_list_i[i_parity + 2*(i_mj + num_mj_n_i*(a + b*n_s))], j, pn2, -phase1*phase2);
        }
      }
    }
  }

  for (int j = 1; j <= n_sds_p_f; j++) {
    if (w_from_p(j, n_s, n_proton, w_shell) > w_max_p_f) {continue;}
    float mj = m_from_p(j, n_s, n_proton, jz_shell);
    if ((mj < mj_min_p_f) || (mj > mj_max_p_f)) {continue;}
    int j_min = j_min_from_p(n_s, n_proton, j);
    
    for (int b = j_min - 1; b < n_s; b++) {
      int phase1;
      int pn1 = a_op(n_s, n_proton, j, b + 1, &phase1, j_min);
      if (pn1 == 0) {continue;}
      for (int a = j_min - 1; a < b; a++) {
        int phase2;
        int pn2 = a_op(n_s, n_proton - 1, pn1, a + 1, &phase2, j_min);
        if (pn2 == 0) {continue;}
        mj = m_from_p(pn2, n_s, n_proton - 2, jz_shell);
        if ((mj < mj_min_p_i) || (mj > mj_max_p_i)) {continue;}
        if (w_from_p(pn2, n_s, n_proton - 2, w_shell) > w_max_p_i) {continue;}
        int i_mj = mj - mj_min_p_i;
        int i_parity = (parity_from_p(pn2, n_s, n_proton - 2, l_shell) + 1)/2;

        if (p2_list_f[i_parity + 2*(i_mj + num_mj_p_i*(b + a*n_s))] == NULL) {
          p2_list_f[i_parity + 2*(i_mj + num_mj_p_i*(b + a*n_s))] = create_sd_node(j, pn2, phase1*phase2, NULL);
        } else {
          sd_append(p2_list_f[i_parity + 2*(i_mj + num_mj_p_i*(b + a*n_s))], j, pn2, phase1*phase2);
        }
        if (p2_list_f[i_parity + 2*(i_mj + num_mj_p_i*(a + b*n_s))] == NULL) {
          p2_list_f[i_parity + 2*(i_mj + num_mj_p_i*(a + b*n_s))] = create_sd_node(j, pn2, -phase1*phase2, NULL);
        } else {
          sd_append(p2_list_f[i_parity + 2*(i_mj + num_mj_p_i*(a + b*n_s))], j, pn2, -phase1*phase2);
        }
      }
    }
  }


  return;
}


void build_two_body_jumps_dmtp2(int n_s, int n_proton, int num_mj_p_i, float mj_min_p_i, float mj_max_p_i, int num_mj_p_f, float mj_min_p_f, float mj_max_p_f, int n_sds_p_f, sd_list** p2_list_f, int n_neutron, int num_mj_n_i, float mj_min_n_i, float mj_max_n_i, int num_mj_n_f, float mj_min_n_f, float mj_max_n_f, int n_sds_n_i, sd_list** n2_list_i, int* jz_shell, int* l_shell) {

  for (int j = 1; j <= n_sds_n_i; j++) {
    float mj = m_from_p(j, n_s, n_neutron, jz_shell);
    if ((mj < mj_min_n_i) || (mj > mj_max_n_i)) {continue;}
    int i_mj = mj - mj_min_n_i;
    int i_parity = (parity_from_p(j, n_s, n_neutron, l_shell) + 1)/2;
    int j_min = j_min_from_p(n_s, n_neutron, j);
    
    for (int b = j_min - 1; b < n_s; b++) {
      int phase1;
      int pn1 = a_op(n_s, n_neutron, j, b + 1, &phase1, j_min);
      if (pn1 == 0) {continue;}
      for (int a = j_min - 1; a < b; a++) {
        int phase2;
        int pn2 = a_op(n_s, n_neutron - 1, pn1, a + 1, &phase2, j_min);
        if (pn2 == 0) {continue;}
        mj = m_from_p(pn2, n_s, n_neutron - 2, jz_shell);
        if ((mj < mj_min_n_f) || (mj > mj_max_n_f)) {continue;}
 
        if (n2_list_i[i_parity + 2*(i_mj + num_mj_n_i*(b + a*n_s))] == NULL) {
          n2_list_i[i_parity + 2*(i_mj + num_mj_n_i*(b + a*n_s))] = create_sd_node(j, pn2, phase1*phase2, NULL);
        } else {
          sd_append(n2_list_i[i_parity + 2*(i_mj + num_mj_n_i*(b + a*n_s))], j, pn2, phase1*phase2);
        }
        if (n2_list_i[i_parity + 2*(i_mj + num_mj_n_i*(a + b*n_s))] == NULL) {
          n2_list_i[i_parity + 2*(i_mj + num_mj_n_i*(a + b*n_s))] = create_sd_node(j, pn2, -phase1*phase2, NULL);
        } else {
          sd_append(n2_list_i[i_parity + 2*(i_mj + num_mj_n_i*(a + b*n_s))], j, pn2, -phase1*phase2);
        }
      }
    }
  }

  for (int j = 1; j <= n_sds_p_f; j++) {
    float mj = m_from_p(j, n_s, n_proton, jz_shell);
    if ((mj < mj_min_p_f) || (mj > mj_max_p_f)) {continue;}

    int j_min = j_min_from_p(n_s, n_proton, j);
    for (int b = j_min - 1; b < n_s; b++) {
      int phase1;
      int pn1 = a_op(n_s, n_proton, j, b + 1, &phase1, j_min);
      if (pn1 == 0) {continue;}
      for (int a = j_min - 1; a < b; a++) {
        int phase2;
        int pn2 = a_op(n_s, n_proton - 1, pn1, a + 1, &phase2, j_min);
        if (pn2 == 0) {continue;}
        mj = m_from_p(pn2, n_s, n_proton - 2, jz_shell);
        if ((mj < mj_min_p_i) || (mj > mj_max_p_i)) {continue;}
        int i_mj = mj - mj_min_p_i;
        int i_parity = (parity_from_p(pn2, n_s, n_proton - 2, l_shell) + 1)/2;

        if (p2_list_f[i_parity + 2*(i_mj + num_mj_p_i*(b + a*n_s))] == NULL) {
          p2_list_f[i_parity + 2*(i_mj + num_mj_p_i*(b + a*n_s))] = create_sd_node(j, pn2, phase1*phase2, NULL);
        } else {
          sd_append(p2_list_f[i_parity + 2*(i_mj + num_mj_p_i*(b + a*n_s))], j, pn2, phase1*phase2);
        }
        if (p2_list_f[i_parity + 2*(i_mj + num_mj_p_i*(a + b*n_s))] == NULL) {
          p2_list_f[i_parity + 2*(i_mj + num_mj_p_i*(a + b*n_s))] = create_sd_node(j, pn2, -phase1*phase2, NULL);
        } else {
          sd_append(p2_list_f[i_parity + 2*(i_mj + num_mj_p_i*(a + b*n_s))], j, pn2, -phase1*phase2);
        }
      }
    }
  }


  return;
}

void build_two_body_jumps_dmtp1_trunc(int n_s, int n_proton_i, int num_mj_p_i, float mj_min_p_i, float mj_max_p_i, int n_proton_f, int num_mj_p_f, float mj_min_p_f, float mj_max_p_f, sd_list** p1_list_i, sd_list** p1_list_f, int* p2_array_f, int n_neutron_i, int num_mj_n_i, float mj_min_n_i, float mj_max_n_i, int n_neutron_f, int num_mj_n_f, float mj_min_n_f, float mj_max_n_f, sd_list** n1_list_i, sd_list** n2_list_i, int* n1_array_f, int* jz_shell, int* l_shell, int* w_shell, int w_max_p_i, int w_max_p_f, int w_max_n_i, int w_max_n_f) {

  unsigned int n_sds_n_i = get_num_sds(n_s, n_neutron_i);
  unsigned int n_sds_n_f = get_num_sds(n_s, n_neutron_f);
  unsigned int n_sds_p_i = get_num_sds(n_s, n_proton_i);
  unsigned int n_sds_p_f = get_num_sds(n_s, n_proton_f);

  unsigned int n_sds_n_int1f = get_num_sds(n_s, n_neutron_f - 1);
  unsigned int n_sds_p_int2f = get_num_sds(n_s, n_proton_f - 2);


  for (int j = 1; j <= n_sds_n_i; j++) {
    if (w_from_p(j, n_s, n_neutron_i, w_shell) > w_max_n_i) {continue;}
    float mj = m_from_p(j, n_s, n_neutron_i, jz_shell);
    if ((mj < mj_min_n_i) || (mj > mj_max_n_i)) {continue;}
    int i_mj = mj - mj_min_n_i;
    int i_parity = (parity_from_p(j, n_s, n_neutron_i, l_shell) + 1)/2;
    int j_min = j_min_from_p(n_s, n_neutron_i, j);

    for (int b = j_min - 1; b < n_s; b++) {
      int phase1;
      int pn1 = a_op(n_s, n_neutron_i, j, b + 1, &phase1, j_min);
      if (pn1 == 0) {continue;}
      if (w_from_p(pn1, n_s, n_neutron_f, w_shell) <= w_max_n_f) {
        if (n1_list_i[i_parity + 2*(i_mj + num_mj_n_i*b)] == NULL) {
          n1_list_i[i_parity + 2*(i_mj + num_mj_n_i*b)] = create_sd_node(j, pn1, phase1, NULL);
        } else {
          sd_append(n1_list_i[i_parity + 2*(i_mj + num_mj_n_i*b)], j, pn1, phase1);
        }
      }
      for (int a = j_min - 1; a < b; a++) {
        int phase2;
        int pn2 = a_op(n_s, n_neutron_i - 1, pn1, a + 1, &phase2, j_min);
        if (pn2 == 0) {continue;}
 
        if (n2_list_i[i_parity + 2*(i_mj + num_mj_n_i*(b + a*n_s))] == NULL) {
          n2_list_i[i_parity + 2*(i_mj + num_mj_n_i*(b + a*n_s))] = create_sd_node(j, pn2, phase1*phase2, NULL);
        } else {
          sd_append(n2_list_i[i_parity + 2*(i_mj + num_mj_n_i*(b + a*n_s))], j, pn2, phase1*phase2);
        }
        if (n2_list_i[i_parity + 2*(i_mj + num_mj_n_i*(a + b*n_s))] == NULL) {
          n2_list_i[i_parity + 2*(i_mj + num_mj_n_i*(a + b*n_s))] = create_sd_node(j, pn2, -phase1*phase2, NULL);
        } else {
          sd_append(n2_list_i[i_parity + 2*(i_mj + num_mj_n_i*(a + b*n_s))], j, pn2, -phase1*phase2);
        }
      }
    }
  }

  for (int j = 1; j <= n_sds_n_f; j++) {
    float mj = m_from_p(j, n_s, n_neutron_f, jz_shell);
    if ((mj < mj_min_n_f) || (mj > mj_max_n_f)) {continue;}
    int j_min = j_min_from_p(n_s, n_neutron_f, j);

    for (int b = j_min - 1; b < n_s; b++) {
      int phase1;
      int pn1 = a_op(n_s, n_neutron_f, j, b + 1, &phase1, j_min);
      if (pn1 == 0) {continue;}
      n1_array_f[(pn1 - 1) + b*n_sds_n_int1f] = j*phase1;
    }
  }

  for (int j = 1; j <= n_sds_p_i; j++) {
    if (w_from_p(j, n_s, n_proton_i, w_shell) > w_max_n_i) {continue;}
    float mj = m_from_p(j, n_s, n_proton_i, jz_shell);
    if ((mj < mj_min_p_i) || (mj > mj_max_p_i)) {continue;}
    int i_mj = mj - mj_min_p_i;
    int i_parity = (parity_from_p(j, n_s, n_proton_i, l_shell) + 1)/2;
    int j_min = j_min_from_p(n_s, n_proton_i, j);
    
    for (int b = j_min - 1; b < n_s; b++) {
      int phase1;
      int pn1 = a_op(n_s, n_proton_i, j, b + 1, &phase1, j_min);
      if (pn1 == 0) {continue;}
      if (p1_list_i[i_parity + 2*(i_mj + num_mj_p_i*b)] == NULL) {
        p1_list_i[i_parity + 2*(i_mj + num_mj_p_i*b)] = create_sd_node(j, pn1, phase1, NULL);
      } else {
        sd_append(p1_list_i[i_parity + 2*(i_mj + num_mj_p_i*b)], j, pn1, phase1);
      }
    }
  }


  for (int j = 1; j <= n_sds_p_f; j++) {
    if (w_from_p(j, n_s, n_proton_f, w_shell) > w_max_p_f) {continue;}
    float mj = m_from_p(j, n_s, n_proton_f, jz_shell);
    if ((mj < mj_min_p_f) || (mj > mj_max_p_f)) {continue;}
    int j_min = j_min_from_p(n_s, n_proton_f, j);
    
    for (int b = j_min - 1; b < n_s; b++) {
      int phase1;
      int pn1 = a_op(n_s, n_proton_f, j, b + 1, &phase1, j_min);
      if (pn1 == 0) {continue;}
      
      if (w_from_p(pn1, n_s, n_proton_i, w_shell) <= w_max_p_f) {
        mj = m_from_p(pn1, n_s, n_proton_i, jz_shell);
        if ((mj >= mj_min_p_i) && (mj <= mj_max_p_i)) {
          int i_mj = mj - mj_min_p_i;
          int i_parity = (parity_from_p(pn1, n_s, n_proton_i, l_shell) + 1)/2;

          if (p1_list_f[i_parity + 2*(i_mj + num_mj_p_i*b)] == NULL) {
            p1_list_f[i_parity + 2*(i_mj + num_mj_p_i*b)] = create_sd_node(pn1, j, phase1, NULL);
          } else {
            sd_append(p1_list_f[i_parity + 2*(i_mj + num_mj_p_i*b)], pn1, j, phase1);
          }
	}
      }
      for (int a = j_min - 1; a < b; a++) {
        int phase2;
        int pn2 = a_op(n_s, n_proton_f - 1, pn1, a + 1, &phase2, j_min);
        if (pn2 == 0) {continue;}
        p2_array_f[(pn2 - 1) + n_sds_p_int2f*(b + a*n_s)] = phase1*phase2*j;
        p2_array_f[(pn2 - 1) + n_sds_p_int2f*(a + b*n_s)] = -phase1*phase2*j;
      }
    }
  }

  return;
}


void build_two_body_jumps_dmtp1(int n_s, int n_proton_i, int num_mj_p_i, float mj_min_p_i, float mj_max_p_i, int n_proton_f, int num_mj_p_f, float mj_min_p_f, float mj_max_p_f, sd_list** p1_list_i, sd_list** p1_list_f, int* p2_array_f, int n_neutron_i, int num_mj_n_i, float mj_min_n_i, float mj_max_n_i, int n_neutron_f, int num_mj_n_f, float mj_min_n_f, float mj_max_n_f, sd_list** n1_list_i, sd_list** n2_list_i, int* n1_array_f, int* jz_shell, int* l_shell) {

  unsigned int n_sds_n_i = get_num_sds(n_s, n_neutron_i);
  unsigned int n_sds_n_f = get_num_sds(n_s, n_neutron_f);
  unsigned int n_sds_p_i = get_num_sds(n_s, n_proton_i);
  unsigned int n_sds_p_f = get_num_sds(n_s, n_proton_f);

  unsigned int n_sds_n_int1f = get_num_sds(n_s, n_neutron_f - 1);
  unsigned int n_sds_p_int2f = get_num_sds(n_s, n_proton_f - 2);


  for (int j = 1; j <= n_sds_n_i; j++) {
    float mj = m_from_p(j, n_s, n_neutron_i, jz_shell);
    if ((mj < mj_min_n_i) || (mj > mj_max_n_i)) {continue;}
    int i_mj = mj - mj_min_n_i;
    int i_parity = (parity_from_p(j, n_s, n_neutron_i, l_shell) + 1)/2;
    int j_min = j_min_from_p(n_s, n_neutron_i, j);
    
    for (int b = j_min - 1; b < n_s; b++) {
      int phase1;
      int pn1 = a_op(n_s, n_neutron_i, j, b + 1, &phase1, j_min);
      if (pn1 == 0) {continue;}
      if (n1_list_i[i_parity + 2*(i_mj + num_mj_n_i*b)] == NULL) {
        n1_list_i[i_parity + 2*(i_mj + num_mj_n_i*b)] = create_sd_node(j, pn1, phase1, NULL);
      } else {
        sd_append(n1_list_i[i_parity + 2*(i_mj + num_mj_n_i*b)], j, pn1, phase1);
      }

      for (int a = j_min - 1; a < b; a++) {
        int phase2;
        int pn2 = a_op(n_s, n_neutron_i - 1, pn1, a + 1, &phase2, j_min);
        if (pn2 == 0) {continue;}
 
        if (n2_list_i[i_parity + 2*(i_mj + num_mj_n_i*(b + a*n_s))] == NULL) {
          n2_list_i[i_parity + 2*(i_mj + num_mj_n_i*(b + a*n_s))] = create_sd_node(j, pn2, phase1*phase2, NULL);
        } else {
          sd_append(n2_list_i[i_parity + 2*(i_mj + num_mj_n_i*(b + a*n_s))], j, pn2, phase1*phase2);
        }
        if (n2_list_i[i_parity + 2*(i_mj + num_mj_n_i*(a + b*n_s))] == NULL) {
          n2_list_i[i_parity + 2*(i_mj + num_mj_n_i*(a + b*n_s))] = create_sd_node(j, pn2, -phase1*phase2, NULL);
        } else {
          sd_append(n2_list_i[i_parity + 2*(i_mj + num_mj_n_i*(a + b*n_s))], j, pn2, -phase1*phase2);
        }
      }
    }
  }

  for (int j = 1; j <= n_sds_n_f; j++) {
    float mj = m_from_p(j, n_s, n_neutron_f, jz_shell);
    if ((mj < mj_min_n_f) || (mj > mj_max_n_f)) {continue;}
    int j_min = j_min_from_p(n_s, n_neutron_f, j);
    for (int b = j_min - 1; b < n_s; b++) {
      int phase1;
      int pn1 = a_op(n_s, n_neutron_f, j, b + 1, &phase1, j_min);
      if (pn1 == 0) {continue;}
      n1_array_f[(pn1 - 1) + b*n_sds_n_int1f] = j*phase1;
    }
  }
  for (int j = 1; j <= n_sds_p_i; j++) {
    float mj = m_from_p(j, n_s, n_proton_i, jz_shell);
    if ((mj < mj_min_p_i) || (mj > mj_max_p_i)) {continue;}
    int i_mj = mj - mj_min_p_i;
    int i_parity = (parity_from_p(j, n_s, n_proton_i, l_shell) + 1)/2;
    int j_min = j_min_from_p(n_s, n_proton_i, j);

    for (int b = j_min - 1; b < n_s; b++) {
      int phase1;
      int pn1 = a_op(n_s, n_proton_i, j, b + 1, &phase1, j_min);
      if (pn1 == 0) {continue;}
      if (p1_list_i[i_parity + 2*(i_mj + num_mj_p_i*b)] == NULL) {
        p1_list_i[i_parity + 2*(i_mj + num_mj_p_i*b)] = create_sd_node(j, pn1, phase1, NULL);
      } else {
        sd_append(p1_list_i[i_parity + 2*(i_mj + num_mj_p_i*b)], j, pn1, phase1);
      }
    }
  }


  for (int j = 1; j <= n_sds_p_f; j++) {
    float mj = m_from_p(j, n_s, n_proton_f, jz_shell);
    if ((mj < mj_min_p_f) || (mj > mj_max_p_f)) {continue;}
    int j_min = j_min_from_p(n_s, n_proton_f, j);

    for (int b = j_min - 1; b < n_s; b++) {
      int phase1;
      int pn1 = a_op(n_s, n_proton_f, j, b + 1, &phase1, j_min);
      if (pn1 == 0) {continue;}
      
      mj = m_from_p(pn1, n_s, n_proton_i, jz_shell);
      if ((mj >= mj_min_p_i) && (mj <= mj_max_p_i)) {
        int i_mj = mj - mj_min_p_i;
        int i_parity = (parity_from_p(pn1, n_s, n_proton_i, l_shell) + 1)/2;

        if (p1_list_f[i_parity + 2*(i_mj + num_mj_p_i*b)] == NULL) {
          p1_list_f[i_parity + 2*(i_mj + num_mj_p_i*b)] = create_sd_node(pn1, j, phase1, NULL);
        } else {
          sd_append(p1_list_f[i_parity + 2*(i_mj + num_mj_p_i*b)], pn1, j, phase1);
        }
      }
      for (int a = j_min - 1; a < b; a++) {
        int phase2;
        int pn2 = a_op(n_s, n_proton_f - 1, pn1, a + 1, &phase2, j_min);
        if (pn2 == 0) {continue;}
        p2_array_f[(pn2 - 1) + n_sds_p_int2f*(b + a*n_s)] = phase1*phase2*j;
        p2_array_f[(pn2 - 1) + n_sds_p_int2f*(a + b*n_s)] = -phase1*phase2*j;
      }
    }
  }

  return;
}

void build_one_body_jumps_dmt0(int n_s, int n_p, float mj_min, float mj_max, int num_mj, int n_sds_i, int n_sds_int, int*a1_array_f, wf_list** a0_list_i, sd_list** a1_list_i, int* jz_shell, int* l_shell) {
/*
  Input(s):
    mj_min_i: 
*/
  for (int j = 1; j <= n_sds_i; j++) {
    float mj = m_from_p(j, n_s, n_p, jz_shell);
    if ((mj < mj_min) || (mj > mj_max)) {continue;}
    int i_parity = (parity_from_p(j, n_s, n_p, l_shell) + 1)/2;
    int i_mj = mj - mj_min;
    if (a0_list_i[i_parity + 2*i_mj] == NULL) {
      a0_list_i[i_parity + 2*i_mj] = create_wf_node(j, NULL);
    } else {
      wf_append(a0_list_i[i_parity + 2*i_mj], j);
    }
    int j_min = j_min_from_p(n_s, n_p, j);
    
    for (int a = j_min - 1; a < n_s; a++) {
      int phase;
      int pn = a_op(n_s, n_p, j, a + 1, &phase, j_min);
      if (pn == 0) {continue;}
      a1_array_f[(pn - 1) + a*n_sds_int] = j*phase;
      if (a1_list_i[i_parity + 2*(i_mj + num_mj*a)] == NULL) {
        a1_list_i[i_parity + 2*(i_mj + num_mj*a)] = create_sd_node(j, pn, phase, NULL);
      } else {
        sd_append(a1_list_i[i_parity + 2*(i_mj + num_mj*a)], j, pn, phase);
      }
    }
  } 

  return;
}

void build_one_body_jumps_dmtp1(int n_s, int n_proton_f, int num_mj_p_i, float mj_min_p_i, float mj_max_p_i, int num_mj_p_f, float mj_min_p_f, float mj_max_p_f, sd_list** p1_list_f, int n_neutron_i, int num_mj_n_i, float mj_min_n_i, float mj_max_n_i, sd_list** n1_list_i, int* jz_shell, int* l_shell) {

  unsigned int n_sds_n_i = get_num_sds(n_s, n_neutron_i);
  unsigned int n_sds_p_f = get_num_sds(n_s, n_proton_f);

  for (int j = 1; j <= n_sds_n_i; j++) {
    float mj = m_from_p(j, n_s, n_neutron_i, jz_shell);
    if ((mj < mj_min_n_i) || (mj > mj_max_n_i)) {continue;}
    int i_mj = mj - mj_min_n_i;
    int i_parity = (parity_from_p(j, n_s, n_neutron_i, l_shell) + 1)/2;
    int j_min = j_min_from_p(n_s, n_neutron_i, j);

    for (int a = j_min - 1; a < n_s; a++) {
      int phase;
      int pn = a_op(n_s, n_neutron_i, j, a + 1, &phase, j_min);
      if (pn == 0) {continue;}

      if (n1_list_i[i_parity + 2*(i_mj + num_mj_n_i*a)] == NULL) {
        n1_list_i[i_parity + 2*(i_mj + num_mj_n_i*a)] = create_sd_node(j, pn, phase, NULL);
      } else {
        sd_append(n1_list_i[i_parity + 2*(i_mj + num_mj_n_i*a)], j, pn, phase);
      }
    }
  }

  for (int j = 1; j <= n_sds_p_f; j++) {
    float mj = m_from_p(j, n_s, n_proton_f, jz_shell);
    if ((mj < mj_min_p_f) || (mj > mj_max_p_f)) {continue;}
    int j_min = j_min_from_p(n_s, n_proton_f, j);

    for (int a = j_min - 1; a < n_s; a++) {
      int phase;
      int pn = a_op(n_s, n_proton_f, j, a + 1, &phase, j_min);
      if (pn == 0) {continue;}
      mj = m_from_p(pn, n_s, n_proton_f - 1, jz_shell);
      if ((mj < mj_min_p_i) || (mj > mj_max_p_i)) {continue;}
      int i_mj = mj - mj_min_p_i;
      int i_parity = (parity_from_p(pn, n_s, n_proton_f - 1, l_shell) + 1)/2;

      if (p1_list_f[i_parity + 2*(i_mj + num_mj_p_i*a)] == NULL) {
        p1_list_f[i_parity + 2*(i_mj + num_mj_p_i*a)] = create_sd_node(pn, j, phase, NULL);
      } else {
        sd_append(p1_list_f[i_parity + 2*(i_mj + num_mj_p_i*a)], pn, j, phase);
      }
    }
  }

  return;
}



void build_one_body_jumps_dmt0_trunc(int n_s, int n_p, float mj_min, float mj_max, int num_mj, int n_sds_i, int n_sds_int, int*a1_array_f, wf_list** a0_list_i, sd_list** a1_list_i, int* jz_shell, int* l_shell, int* w_shell, int w_max) {
/*
  Input(s):
    mj_min_i: 
*/
  for (int j = 1; j <= n_sds_i; j++) {
    if (w_from_p(j, n_s, n_p, w_shell) > w_max) {continue;}
    float mj = m_from_p(j, n_s, n_p, jz_shell);
    if ((mj < mj_min) || (mj > mj_max)) {continue;}
    int i_parity = (parity_from_p(j, n_s, n_p, l_shell) + 1)/2;
    int i_mj = mj - mj_min;
    if (a0_list_i[i_parity + 2*i_mj] == NULL) {
      a0_list_i[i_parity + 2*i_mj] = create_wf_node(j, NULL);
    } else {
      wf_append(a0_list_i[i_parity + 2*i_mj], j);
    }
    int j_min = j_min_from_p(n_s, n_p, j);
    
    for (int a = j_min - 1; a < n_s; a++) {
      int phase;
      int pn = a_op(n_s, n_p, j, a + 1, &phase, j_min);
      if (pn == 0) {continue;}
      a1_array_f[(pn - 1) + a*n_sds_int] = j*phase;
      if (a1_list_i[i_parity + 2*(i_mj + num_mj*a)] == NULL) {
        a1_list_i[i_parity + 2*(i_mj + num_mj*a)] = create_sd_node(j, pn, phase, NULL);
      } else {
        sd_append(a1_list_i[i_parity + 2*(i_mj + num_mj*a)], j, pn, phase);
      }
    }
  } 

  return;
}

void build_one_body_jumps_dmtp1_trunc(int n_s, int n_proton_f, int num_mj_p_i, float mj_min_p_i, float mj_max_p_i, int num_mj_p_f, float mj_min_p_f, float mj_max_p_f, sd_list** p1_list_f, int n_neutron_i, int num_mj_n_i, float mj_min_n_i, float mj_max_n_i, sd_list** n1_list_i, int* jz_shell, int* l_shell, int* w_shell, int w_max_p_i, int w_max_p_f, int w_max_n_i, int w_max_n_f) {

  unsigned int n_sds_n_i = get_num_sds(n_s, n_neutron_i);
  unsigned int n_sds_p_f = get_num_sds(n_s, n_proton_f);

  for (int j = 1; j <= n_sds_n_i; j++) {
    if (w_from_p(j, n_s, n_neutron_i, w_shell) > w_max_n_i) {continue;}
    float mj = m_from_p(j, n_s, n_neutron_i, jz_shell);
    if ((mj < mj_min_n_i) || (mj > mj_max_n_i)) {continue;}
    int i_mj = mj - mj_min_n_i;
    int i_parity = (parity_from_p(j, n_s, n_neutron_i, l_shell) + 1)/2;
    int j_min = j_min_from_p(n_s, n_neutron_i, j);

    for (int a = j_min - 1; a < n_s; a++) {
      int phase;
      int pn = a_op(n_s, n_neutron_i, j, a + 1, &phase, j_min);
      if (pn == 0) {continue;}
      if (w_from_p(pn, n_s, n_neutron_i - 1, w_shell) > w_max_n_f) {continue;}
      if (n1_list_i[i_parity + 2*(i_mj + num_mj_n_i*a)] == NULL) {
        n1_list_i[i_parity + 2*(i_mj + num_mj_n_i*a)] = create_sd_node(j, pn, phase, NULL);
      } else {
        sd_append(n1_list_i[i_parity + 2*(i_mj + num_mj_n_i*a)], j, pn, phase);
      }
    }
  }

  for (int j = 1; j <= n_sds_p_f; j++) {
    if (w_from_p(j, n_s, n_proton_f, w_shell) > w_max_p_f) {continue;}
    float mj = m_from_p(j, n_s, n_proton_f, jz_shell);
    if ((mj < mj_min_p_f) || (mj > mj_max_p_f)) {continue;}
    int j_min = j_min_from_p(n_s, n_proton_f, j);

    for (int a = j_min - 1; a < n_s; a++) {
      int phase;
      int pn = a_op(n_s, n_proton_f, j, a + 1, &phase, j_min);
      if (pn == 0) {continue;}
      if (w_from_p(pn, n_s, n_proton_f - 1, w_shell) > w_max_p_i) {continue;}
      mj = m_from_p(pn, n_s, n_proton_f - 1, jz_shell);
      if ((mj < mj_min_p_i) || (mj > mj_max_p_i)) {continue;}
      int i_mj = mj - mj_min_p_i;
      int i_parity = (parity_from_p(pn, n_s, n_proton_f - 1, l_shell) + 1)/2;

      if (p1_list_f[i_parity + 2*(i_mj + num_mj_p_i*a)] == NULL) {
        p1_list_f[i_parity + 2*(i_mj + num_mj_p_i*a)] = create_sd_node(pn, j, phase, NULL);
      } else {
        sd_append(p1_list_f[i_parity + 2*(i_mj + num_mj_p_i*a)], pn, j, phase);
      }
    }
  }


  return;
}

void trace_one_body_nodes_dmt0(int a, int b, int num_mj_p_i, float mj_min_p_i, int num_mj_n_i, float mj_min_n_i, int n_sds_int, int* p1_array_f, sd_list** p1_list_i, wf_list** n0_list_i, wfnData* wd, int i_op, eigen_list* transition, double* density) {

  for (int imjp = 0; imjp < num_mj_p_i; imjp++) {
    float mjp = imjp + mj_min_p_i;
    for (int imjn = 0; imjn < num_mj_n_i; imjn++) {
      float mjn = imjn + mj_min_n_i;
      if ((mjp + mjn != 0) && (mjp + mjn != 0.5)) {continue;}
      for (int ipar_p = 0; ipar_p <= 1; ipar_p++) {
        int ipar_n;
	if (wd->parity_i == '+') {
	  if (ipar_p == 0) {
		  ipar_n = 0;
	  } else {
		  ipar_n = 1;
	  }
	} else if (wd->parity_i == '-') {
		if (ipar_p == 0) {
			ipar_n = 1;
		} else {
			ipar_n = 0;
		}
	} else {printf("Parity error\n"); exit(0);}

        sd_list* node1 = p1_list_i[ipar_p + 2*(imjp + num_mj_p_i*b)];
        while (node1 != NULL) {
          unsigned int ppn = node1->pn;
          unsigned int ppi = node1->pi;
          int phase1 = node1->phase;
          int ppf = p1_array_f[(ppn - 1) + n_sds_int*a];
          if (ppf == 0) {node1 = node1->next; continue;}
          int phase2 = 1;
          node1 = node1->next;
          if (ppf < 0) {
            ppf *= -1;
	    phase2 = -1;
          }

          wf_list* node2 = n0_list_i[ipar_n + 2*imjn];
          while (node2 != NULL) {
	    int pn = node2->p;
   	    int index_i = -1;
	    int index_f = -1;

  	    if (i_op == 0) {
	      unsigned int p_hash_i = ppi + wd->n_sds_p_i*(pn % HASH_SIZE);
	      wh_list* node3 = wd->wh_hash_i[p_hash_i];
	      while (node3 != NULL) {
	        if ((pn == node3->pn) && (ppi == node3->pp)) {
	          index_i = node3->index;
	          break;
	        }
	        node3 = node3->next;
	      }
	      if (index_i < 0) {node2 = node2->next; continue;}

	      unsigned int p_hash_f = ppf + wd->n_sds_p_f*(pn % HASH_SIZE);
	      node3 = wd->wh_hash_f[p_hash_f];
	      while (node3 != NULL) {
	        if ((pn == node3->pn) && (ppf == node3->pp)) {
	          index_f = node3->index;
	          break;
	        }
	        node3 = node3->next;
	      }
	      if (index_f < 0) {node2 = node2->next; continue;}
              eigen_list* eig_pair = transition;
              int i_trans = 0;
              while (eig_pair != NULL) {
                int psi_i = eig_pair->eig_i;
                int psi_f = eig_pair->eig_f;
	        density[i_trans] += wd->bc_i[psi_i + wd->n_eig_i*index_i]*wd->bc_f[psi_f + wd->n_eig_f*index_f]*phase1*phase2;
                i_trans++;
                eig_pair = eig_pair->next;
              }
	    } else {
	      unsigned int p_hash_i = pn + wd->n_sds_p_i*(ppi % HASH_SIZE);
	      wh_list* node3 = wd->wh_hash_i[p_hash_i];
	      while (node3 != NULL) {
	        if ((pn == node3->pp) && (ppi == node3->pn)) {
	          index_i = node3->index;
	          break;
	        }
	        node3 = node3->next;
	      }
	      if (index_i < 0) {node2 = node2->next; continue;}

	      unsigned int p_hash_f = pn + wd->n_sds_p_f*(ppf % HASH_SIZE);
	      node3 = wd->wh_hash_f[p_hash_f];
	      while (node3 != NULL) {
	        if ((pn == node3->pp) && (ppf == node3->pn)) {
	          index_f = node3->index;
	          break;
	        }
	        node3 = node3->next;
	      }
	      if (index_f < 0) {node2 = node2->next; continue;}
              eigen_list* eig_pair = transition;
              int i_trans = 0;
              while (eig_pair != NULL) {
                int psi_i = eig_pair->eig_i;
                int psi_f = eig_pair->eig_f;
	        density[i_trans] += wd->bc_i[psi_i + wd->n_eig_i*index_i]*wd->bc_f[psi_f + wd->n_eig_f*index_f]*phase1*phase2;
                i_trans++;
                eig_pair = eig_pair->next;
              }
	    } 
	    node2 = node2->next;
	  }
        } 
      }
    }
  }
  return;
}

void trace_one_body_nodes_dmtp1(int a, int b, int num_mj_p_i, float mj_min_p_i, int num_mj_n_i, float mj_min_n_i, sd_list** n1_list_i, sd_list** p1_list_f, wfnData* wd, int i_op, eigen_list *transition, double* density) {
/* 

*/
  for (int imjp = 0; imjp < num_mj_p_i; imjp++) {
    float mjp = imjp + mj_min_p_i;
    for (int imjn = 0; imjn < num_mj_n_i; imjn++) {
      float mjn = imjn + mj_min_n_i;
      if ((mjp + mjn != 0) && (mjp + mjn != 0.5)) {continue;}
      for (int ipar_p = 0; ipar_p <= 1; ipar_p++) {
        int ipar_n;
	if (wd->parity_i == '+') {
	  if (ipar_p == 0) {
		  ipar_n = 0;
	  } else {
		  ipar_n = 1;
	  }
	} else if (wd->parity_i == '-') {
		if (ipar_p == 0) {
			ipar_n = 1;
		} else {
			ipar_n = 0;
		}
	} else {printf("Parity error\n"); exit(0);}
        sd_list* node1 = n1_list_i[ipar_n + 2*(imjn + num_mj_n_i*b)];
        // Loop over final states resulting from 2x a_op
        while (node1 != NULL) {
          unsigned int pnf = node1->pn; // Get final state p_f
          unsigned int pni = node1->pi; // Get initial state p_i
          int phase1 = node1->phase;
          sd_list* node2 = p1_list_f[ipar_p + 2*(imjp + num_mj_p_i*a)]; //hash corresponds to a_op operators

        // Loop over n_f  
        //int m_pf = m_from_p(ppf, ns, npp, wd->jz_shell);
        while (node2 != NULL) {
          unsigned int ppi = node2->pi;
          unsigned int ppf = node2->pn;
          int phase2 = node2->phase;
          int index_i = -1;
          int index_f = -1;  
          if (i_op == 0) {
	    unsigned int p_hash_i = ppi + wd->n_sds_p_i*(pni % HASH_SIZE);
	    wh_list* node3 = wd->wh_hash_i[p_hash_i];
	    while (node3 != NULL) {
	      if ((pni == node3->pn) && (ppi == node3->pp)) {
	        index_i = node3->index;
	        break;
	      }
	      node3 = node3->next;
	    }
	    if (index_i < 0) {node2 = node2->next; continue;}

	    unsigned int p_hash_f = ppf + wd->n_sds_p_f*(pnf % HASH_SIZE);
	    node3 = wd->wh_hash_f[p_hash_f];
	    while (node3 != NULL) {
	      if ((pnf == node3->pn) && (ppf == node3->pp)) {
	        index_f = node3->index;
	        break;
	      }
	      node3 = node3->next;
	    }
	    if (index_f < 0) {node2 = node2->next; continue;}
            eigen_list* eig_pair = transition;
            int i_trans = 0;
            while (eig_pair != NULL) {
              int psi_i = eig_pair->eig_i;
              int psi_f = eig_pair->eig_f;
	      density[i_trans] += wd->bc_i[psi_i + wd->n_eig_i*index_i]*wd->bc_f[psi_f + wd->n_eig_f*index_f]*phase1*phase2;
              i_trans++;
              eig_pair = eig_pair->next;
            }
	  } else {
	    unsigned int p_hash_i = pni + wd->n_sds_p_i*(ppi % HASH_SIZE);
	    wh_list* node3 = wd->wh_hash_i[p_hash_i];
	    while (node3 != NULL) {
	      if ((pni == node3->pp) && (ppi == node3->pn)) {
	        index_i = node3->index;
	        break;
	      }
	      node3 = node3->next;
	    }
	    if (index_i < 0) {node2 = node2->next; continue;}

	    unsigned int p_hash_f = pnf + wd->n_sds_p_f*(ppf % HASH_SIZE);
	    node3 = wd->wh_hash_f[p_hash_f];
	    while (node3 != NULL) {
	      if ((pnf == node3->pp) && (ppf == node3->pn)) {
	        index_f = node3->index;
	        break;
	      }
	      node3 = node3->next;
	    }
	    if (index_f < 0) {node2 = node2->next; continue;}
            eigen_list* eig_pair = transition;
            int i_trans = 0;
            while (eig_pair != NULL) {
              int psi_i = eig_pair->eig_i;
              int psi_f = eig_pair->eig_f;
	      density[i_trans] += wd->bc_i[psi_i + wd->n_eig_i*index_i]*wd->bc_f[psi_f + wd->n_eig_f*index_f]*phase1*phase2;
              i_trans++;
              eig_pair = eig_pair->next;
            }
	  } 
          node2 = node2->next;
	}
        node1 = node1->next;
	}
      }
    }  
  }
  return;
}

int test_suite() {
  int pass = 1;
  if (!test_log_factorial()) {pass = 0;}
  if (!test_clebsch_gordan()) {pass = 0;}
  if (!test_n_choose_k()) {pass = 0;}
  if (!test_slater_routines()) {pass = 0;}

  if (pass == 0) {printf("Error. One of more unit test failed.\n"); exit(0);}
  printf("\n");
  speedParams* sp = read_parameter_file("../examples/param_files/ne20_ne20_1body_test.param");
  one_body_density(sp);

  FILE *in_file;
  in_file = fopen("../examples/output/ne20_ne20_1body_0_0.dens", "r");
  int Ji, Ti, Jf, Tf, Jop, Top;
  fscanf(in_file, "%d %d %d %d %d %d\n", &Jf, &Tf, &Ji, &Ti, &Jop, &Top);
  if (Jop != 0 || Top != 0) {printf("Error in density matrix\n"); exit(0);}
  double mat = 0;
  int na, ja, nb, jb;
  float dm;
  for (int i = 0; i < 3; i++) {
   fscanf(in_file, "%d %d %d %d %f", &na, &ja, &nb, &jb, &dm);
   mat += sqrt(2.0)*sqrt(ja + 1.0)*dm;
  }
  mat *= 1.0/sqrt((Ji + 1.0)*(Ti + 1.0));
  if (fabs(mat - 4) > pow(10, -4)) {pass = 0;} else {
  printf("Ne20 state #1 -> Ne20 state #1: 1-body J=0 T=0 sum rule: Pass\n");}

  fclose(in_file);


  in_file = fopen("../examples/output/ne20_ne20_1body_1_1.dens", "r");
  fscanf(in_file, "%d %d %d %d %d %d\n", &Jf, &Tf, &Ji, &Ti, &Jop, &Top);
  if (Jop != 0 || Top != 0) {printf("Error in density matrix\n"); exit(0);}
  mat = 0;
  for (int i = 0; i < 3; i++) {
   fscanf(in_file, "%d %d %d %d %f", &na, &ja, &nb, &jb, &dm);
   mat += sqrt(2.0)*sqrt(ja + 1.0)*dm;
  }
  mat *= 1.0/sqrt((Ji + 1.0)*(Ti + 1.0));
  if (fabs(mat - 4) > pow(10, -4)) {pass = 0;} else {
  printf("Ne20 state #2 -> Ne20 state #2: 1-body J=0 T=0 sum rule: Pass\n");}

  fclose(in_file);

  in_file = fopen("../examples/output/ne20_ne20_1body_2_2.dens", "r");
  fscanf(in_file, "%d %d %d %d %d %d\n", &Jf, &Tf, &Ji, &Ti, &Jop, &Top);
  if (Jop != 0 || Top != 0) {printf("Error in density matrix\n"); exit(0);}
  mat = 0;
  for (int i = 0; i < 3; i++) {
   fscanf(in_file, "%d %d %d %d %f", &na, &ja, &nb, &jb, &dm);
   mat += sqrt(2.0)*sqrt(ja + 1.0)*dm;
  }
  mat *= 1.0/sqrt((Ji + 1.0)*(Ti + 1.0));
  if (fabs(mat - 4) > pow(10, -4)) {pass = 0;} else {
  printf("Ne20 state #3 -> Ne20 state #3: 1-body J=0 T=0 sum rule: Pass\n");}
  printf("\n");
  fclose(in_file);

  sp = read_parameter_file("../examples/param_files/ne20_ne20_2body_test.param");
  two_body_density(sp);
  
  mat = compute_2body_J0_T0_sum_rule("../examples/output/ne20_ne20_2body_J0_T0_0_0.dens");

  if (fabs(mat - 6) < pow(10, -4)) {printf("Ne20 state #1 -> Ne20 state #1: 2-body J=0 T=0 Sum Rule: Pass\n");}

  mat = compute_2body_J0_T0_sum_rule("../examples/output/ne20_ne20_2body_J0_T0_1_1.dens")/sqrt(5.0);

  if (fabs(mat - 6) < pow(10, -4)) {printf("Ne20 state #2 -> Ne20 state #2: 2-body J=0 T=0 Sum Rule: Pass\n");}

  mat = compute_2body_J0_T0_sum_rule("../examples/output/ne20_ne20_2body_J0_T0_2_2.dens")/sqrt(9.0);

  if (fabs(mat - 6) < pow(10, -4)) {printf("Ne20 state #3 -> Ne20 state #3: 2-body J=0 T=0 Sum Rule: Pass\n");}
  printf("\n");
  free(sp);

  // Begin Na20 -> Na20 tests

  sp = read_parameter_file("../examples/param_files/na20_na20_1body_test.param");
  one_body_density(sp);

  in_file = fopen("../examples/output/na20_na20_1body_0_0.dens", "r");
  fscanf(in_file, "%d %d %d %d %d %d\n", &Jf, &Tf, &Ji, &Ti, &Jop, &Top);
  if (Jop != 0 || Top != 0) {printf("Error in density matrix\n"); exit(0);}
  mat = 0;
  for (int i = 0; i < 3; i++) {
   fscanf(in_file, "%d %d %d %d %f", &na, &ja, &nb, &jb, &dm);
   mat += sqrt(2.0)*sqrt(ja + 1.0)*dm;
  }
  mat *= 1.0/sqrt((Ji + 1.0)*(Ti + 1.0));
  if (fabs(mat - 4) > pow(10, -4)) {pass = 0;} else {
  printf("Na20 state #1 -> Na20 state #1: 1-body J=0 T=0 sum rule: Pass\n");}
  fscanf(in_file, "%*d\n");
  fscanf(in_file, "%d %d %d %d %d %d\n", &Jf, &Tf, &Ji, &Ti, &Jop, &Top);
  if (Jop != 0 || Top != 2) {printf("Error in density matrix\n"); exit(0);}
  mat = 0;
  for (int i = 0; i < 3; i++) {
   fscanf(in_file, "%d %d %d %d %f", &na, &ja, &nb, &jb, &dm);
   mat += sqrt(6.0)*sqrt(ja + 1.0)*dm;
  }
  mat *= 1.0/sqrt((Ji + 1.0)*(Ti + 1.0))*(1.0/sqrt(2.0));
  if (fabs(mat - 2) > pow(10, -4)) {pass = 0;} else {
  printf("Na20 state #1 -> Na20 state #1: 1-body J=0 T=1 sum rule: Pass\n");}

  fclose(in_file);

  in_file = fopen("../examples/output/na20_na20_1body_1_1.dens", "r");
  fscanf(in_file, "%d %d %d %d %d %d\n", &Jf, &Tf, &Ji, &Ti, &Jop, &Top);
  if (Jop != 0 || Top != 0) {printf("Error in density matrix\n"); exit(0);}
  mat = 0;
  for (int i = 0; i < 3; i++) {
   fscanf(in_file, "%d %d %d %d %f", &na, &ja, &nb, &jb, &dm);
   mat += sqrt(2.0)*sqrt(ja + 1.0)*dm;
  }
  mat *= 1.0/sqrt((Ji + 1.0)*(Ti + 1.0));
  if (fabs(mat - 4) > pow(10, -4)) {pass = 0;} else {
  printf("Na20 state #2 -> Na20 state #2: 1-body J=0 T=0 sum rule: Pass\n");}
  fscanf(in_file, "%*d\n");
  fscanf(in_file, "%d %d %d %d %d %d\n", &Jf, &Tf, &Ji, &Ti, &Jop, &Top);
  if (Jop != 0 || Top != 2) {printf("Error in density matrix\n"); exit(0);}
  mat = 0;
  for (int i = 0; i < 3; i++) {
   fscanf(in_file, "%d %d %d %d %f", &na, &ja, &nb, &jb, &dm);
   mat += sqrt(6.0)*sqrt(ja + 1.0)*dm;
  }
  mat *= 1.0/sqrt((Ji + 1.0)*(Ti + 1.0))*(1.0/sqrt(2.0));
  if (fabs(mat - 2) > pow(10, -4)) {pass = 0;} else {
  printf("Na20 state #2 -> Na20 state #2: 1-body J=0 T=1 sum rule: Pass\n");}

  fclose(in_file);

  in_file = fopen("../examples/output/na20_na20_1body_2_2.dens", "r");
  fscanf(in_file, "%d %d %d %d %d %d\n", &Jf, &Tf, &Ji, &Ti, &Jop, &Top);
  if (Jop != 0 || Top != 0) {printf("Error in density matrix\n"); exit(0);}
  mat = 0;
  for (int i = 0; i < 3; i++) {
   fscanf(in_file, "%d %d %d %d %f", &na, &ja, &nb, &jb, &dm);
   mat += sqrt(2.0)*sqrt(ja + 1.0)*dm;
  }
  mat *= 1.0/sqrt((Ji + 1.0)*(Ti + 1.0));
  if (fabs(mat - 4) > pow(10, -4)) {pass = 0;} else {
  printf("Na20 state #3 -> Na20 state #3: 1-body J=0 T=0 sum rule: Pass\n");}
  fscanf(in_file, "%*d\n");
  fscanf(in_file, "%d %d %d %d %d %d\n", &Jf, &Tf, &Ji, &Ti, &Jop, &Top);
  if (Jop != 0 || Top != 2) {printf("Error in density matrix\n"); exit(0);}
  mat = 0;
  for (int i = 0; i < 3; i++) {
   fscanf(in_file, "%d %d %d %d %f", &na, &ja, &nb, &jb, &dm);
   mat += sqrt(6.0)*sqrt(ja + 1.0)*dm;
  }
  mat *= 1.0/sqrt((Ji + 1.0)*(Ti + 1.0))*(1.0/sqrt(2.0));
  if (fabs(mat - 2) > pow(10, -4)) {pass = 0;} else {
  printf("Na20 state #3 -> Na20 state #3: 1-body J=0 T=1 sum rule: Pass\n");}

  fclose(in_file);

  printf("\n");

  sp = read_parameter_file("../examples/param_files/na20_na20_2body_test_T0.param");
  two_body_density(sp);
  

  free(sp);

  sp = read_parameter_file("../examples/param_files/na20_na20_2body_test_T1.param");
  two_body_density(sp);
  
  free(sp);

  mat = compute_2body_J0_T0_sum_rule("../examples/output/na20_na20_2body_J0_T0_0_0.dens")/sqrt(15.0);
  if (fabs(mat - 6) < pow(10, -4)) {printf("Na20 state #1 -> Na20 state #1: 2-body J=0 T=0 Sum Rule: Pass\n");}

  mat = compute_2body_J0_T1_sum_rule("../examples/output/na20_na20_2body_J0_T1_0_0.dens")*(1.0/sqrt(30.0));
  if (fabs(mat - 6) < pow(10, -4)) {printf("Na20 state #1 -> Na20 state #1: 2-body J=0 T=1 Sum Rule: Pass\n");}
  
  
  mat = compute_2body_J0_T0_sum_rule("../examples/output/na20_na20_2body_J0_T0_1_1.dens")/sqrt(21.0);
  if (fabs(mat - 6) < pow(10, -4)) {printf("Na20 state #2 -> Na20 state #2: 2-body J=0 T=0 Sum Rule: Pass\n");}

  mat = compute_2body_J0_T1_sum_rule("../examples/output/na20_na20_2body_J0_T1_1_1.dens")*(1.0/sqrt(42.0));
  if (fabs(mat - 6) < pow(10, -4)) {printf("Na20 state #2 -> Na20 state #2: 2-body J=0 T=1 Sum Rule: Pass\n");}

  mat = compute_2body_J0_T0_sum_rule("../examples/output/na20_na20_2body_J0_T0_2_2.dens")/(3.0*sqrt(3.0));
  if (fabs(mat - 6) < pow(10, -4)) {printf("Na20 state #3 -> Na20 state #3: 2-body J=0 T=0 Sum Rule: Pass\n");}

  mat = compute_2body_J0_T1_sum_rule("../examples/output/na20_na20_2body_J0_T1_2_2.dens")*(1.0/(3*sqrt(6.0)));
  if (fabs(mat - 6) < pow(10, -4)) {printf("Na20 state #3 -> Na20 state #3: 2-body J=0 T=1 Sum Rule: Pass\n");}

  printf("\n");


// Begin F20 -> F20 tests

  sp = read_parameter_file("../examples/param_files/f20_f20_1body_test.param");
  one_body_density(sp);

  in_file = fopen("../examples/output/f20_f20_1body_0_0.dens", "r");
  fscanf(in_file, "%d %d %d %d %d %d\n", &Jf, &Tf, &Ji, &Ti, &Jop, &Top);
  if (Jop != 0 || Top != 0) {printf("Error in density matrix\n"); exit(0);}
  mat = 0;
  for (int i = 0; i < 3; i++) {
   fscanf(in_file, "%d %d %d %d %f", &na, &ja, &nb, &jb, &dm);
   mat += sqrt(2.0)*sqrt(ja + 1.0)*dm;
  }
  mat *= 1.0/sqrt((Ji + 1.0)*(Ti + 1.0));
  if (fabs(mat - 4) > pow(10, -4)) {pass = 0;} else {
  printf("F20  state #1 -> F20  state #1: 1-body J=0 T=0 sum rule: Pass\n");}
  fscanf(in_file, "%*d\n");
  fscanf(in_file, "%d %d %d %d %d %d\n", &Jf, &Tf, &Ji, &Ti, &Jop, &Top);
  if (Jop != 0 || Top != 2) {printf("Error in density matrix\n"); exit(0);}
  mat = 0;
  for (int i = 0; i < 3; i++) {
   fscanf(in_file, "%d %d %d %d %f", &na, &ja, &nb, &jb, &dm);
   mat += sqrt(6.0)*sqrt(ja + 1.0)*dm;
  }
  mat *= 1.0/sqrt((Ji + 1.0)*(Ti + 1.0))*(-1.0/sqrt(2.0));
  if (fabs(mat + 2) > pow(10, -4)) {pass = 0;} else {
  printf("F20  state #1 -> F20  state #1: 1-body J=0 T=1 sum rule: Pass\n");}

  fclose(in_file);


  in_file = fopen("../examples/output/f20_f20_1body_1_1.dens", "r");
  fscanf(in_file, "%d %d %d %d %d %d\n", &Jf, &Tf, &Ji, &Ti, &Jop, &Top);
  if (Jop != 0 || Top != 0) {printf("Error in density matrix\n"); exit(0);}
  mat = 0;
  for (int i = 0; i < 3; i++) {
   fscanf(in_file, "%d %d %d %d %f", &na, &ja, &nb, &jb, &dm);
   mat += sqrt(2.0)*sqrt(ja + 1.0)*dm;
  }
  mat *= 1.0/sqrt((Ji + 1.0)*(Ti + 1.0));
  if (fabs(mat - 4) > pow(10, -4)) {pass = 0;} else {
  printf("F20  state #2 -> F20  state #2: 1-body J=0 T=0 sum rule: Pass\n");}
  fscanf(in_file, "%*d\n");
  fscanf(in_file, "%d %d %d %d %d %d\n", &Jf, &Tf, &Ji, &Ti, &Jop, &Top);
  if (Jop != 0 || Top != 2) {printf("Error in density matrix\n"); exit(0);}
  mat = 0;
  for (int i = 0; i < 3; i++) {
   fscanf(in_file, "%d %d %d %d %f", &na, &ja, &nb, &jb, &dm);
   mat += sqrt(6.0)*sqrt(ja + 1.0)*dm;
  }
  mat *= 1.0/sqrt((Ji + 1.0)*(Ti + 1.0))*(-1.0/sqrt(2.0));
  if (fabs(mat + 2) > pow(10, -4)) {pass = 0;} else {
  printf("F20  state #2 -> F20  state #2: 1-body J=0 T=1 sum rule: Pass\n");}

  fclose(in_file);

  in_file = fopen("../examples/output/f20_f20_1body_2_2.dens", "r");
  fscanf(in_file, "%d %d %d %d %d %d\n", &Jf, &Tf, &Ji, &Ti, &Jop, &Top);
  if (Jop != 0 || Top != 0) {printf("Error in density matrix\n"); exit(0);}
  mat = 0;
  for (int i = 0; i < 3; i++) {
   fscanf(in_file, "%d %d %d %d %f", &na, &ja, &nb, &jb, &dm);
   mat += sqrt(2.0)*sqrt(ja + 1.0)*dm;
  }
  mat *= 1.0/sqrt((Ji + 1.0)*(Ti + 1.0));
  if (fabs(mat - 4) > pow(10, -4)) {pass = 0;} else {
  printf("F20  state #3 -> F20  state #3: 1-body J=0 T=0 sum rule: Pass\n");}
  fscanf(in_file, "%*d\n");
  fscanf(in_file, "%d %d %d %d %d %d\n", &Jf, &Tf, &Ji, &Ti, &Jop, &Top);
  if (Jop != 0 || Top != 2) {printf("Error in density matrix\n"); exit(0);}
  mat = 0;
  for (int i = 0; i < 3; i++) {
   fscanf(in_file, "%d %d %d %d %f", &na, &ja, &nb, &jb, &dm);
   mat += sqrt(6.0)*sqrt(ja + 1.0)*dm;
  }
  mat *= 1.0/sqrt((Ji + 1.0)*(Ti + 1.0))*(-1.0/sqrt(2.0));
  if (fabs(mat + 2) > pow(10, -4)) {pass = 0;} else {
  printf("F20  state #3 -> F20  state #3: 1-body J=0 T=1 sum rule: Pass\n");}

  fclose(in_file);

  printf("\n");

  sp = read_parameter_file("../examples/param_files/f20_f20_2body_test_T0.param");
  two_body_density(sp);
  free(sp);

  sp = read_parameter_file("../examples/param_files/f20_f20_2body_test_T1.param");
  two_body_density(sp);
  free(sp);

  mat = compute_2body_J0_T0_sum_rule("../examples/output/f20_f20_2body_J0_T0_0_0.dens")/sqrt(15.0);
  if (fabs(mat - 6) < pow(10, -4)) {printf("F20  state #1 -> F20  state #1: 2-body J=0 T=0 Sum Rule: Pass\n");}

  mat = compute_2body_J0_T1_sum_rule("../examples/output/f20_f20_2body_J0_T1_0_0.dens")*(-1.0/sqrt(30.0));
  if (fabs(mat + 6) < pow(10, -4)) {printf("F20  state #1 -> F20  state #1: 2-body J=0 T=1 Sum Rule: Pass\n");}

  mat = compute_2body_J0_T0_sum_rule("../examples/output/f20_f20_2body_J0_T0_1_1.dens")/sqrt(21.0);
  if (fabs(mat - 6) < pow(10, -4)) {printf("F20  state #2 -> F20  state #2: 2-body J=0 T=0 Sum Rule: Pass\n");}

  mat = compute_2body_J0_T1_sum_rule("../examples/output/f20_f20_2body_J0_T1_1_1.dens")*(-1.0/sqrt(42.0));
  if (fabs(mat + 6) < pow(10, -4)) {printf("F20  state #2 -> F20  state #2: 2-body J=0 T=1 Sum Rule: Pass\n");}

  mat = compute_2body_J0_T0_sum_rule("../examples/output/f20_f20_2body_J0_T0_2_2.dens")/(3.0*sqrt(3.0));
  if (fabs(mat - 6) < pow(10, -4)) {printf("F20  state #3 -> F20  state #3: 2-body J=0 T=0 Sum Rule: Pass\n");}

  mat = compute_2body_J0_T1_sum_rule("../examples/output/f20_f20_2body_J0_T1_2_2.dens")*(-1.0/(3*sqrt(6.0)));
  if (fabs(mat + 6) < pow(10, -4)) {printf("F20  state #3 -> F20  state #3: 2-body J=0 T=1 Sum Rule: Pass\n");}

  printf("\n");

  // Begin Mg20 -> Mg20 tests

  sp = read_parameter_file("../examples/param_files/mg20_mg20_1body_test.param");
  one_body_density(sp);

  in_file = fopen("../examples/output/mg20_mg20_1body_0_0.dens", "r");
  fscanf(in_file, "%d %d %d %d %d %d\n", &Jf, &Tf, &Ji, &Ti, &Jop, &Top);
  if (Jop != 0 || Top != 0) {printf("Error in density matrix\n"); exit(0);}
  mat = 0;
  for (int i = 0; i < 3; i++) {
   fscanf(in_file, "%d %d %d %d %f", &na, &ja, &nb, &jb, &dm);
   mat += sqrt(2.0)*sqrt(ja + 1.0)*dm;
  }
  mat *= 1.0/sqrt((Ji + 1.0)*(Ti + 1.0));
  if (fabs(mat - 4) > pow(10, -4)) {pass = 0;} else {
  printf("Mg20 state #1 -> Mg20 state #1: 1-body J=0 T=0 sum rule: Pass\n");}
  fscanf(in_file, "%*d\n");
  fscanf(in_file, "%d %d %d %d %d %d\n", &Jf, &Tf, &Ji, &Ti, &Jop, &Top);
  if (Jop != 0 || Top != 2) {printf("Error in density matrix\n"); exit(0);}
  mat = 0;
  for (int i = 0; i < 3; i++) {
   fscanf(in_file, "%d %d %d %d %f", &na, &ja, &nb, &jb, &dm);
   mat += sqrt(6.0)*sqrt(ja + 1.0)*dm;
  }
  mat *= 1.0/sqrt((Ji + 1.0)*(Ti + 1.0))*(sqrt(2.0/3.0));
  if (fabs(mat - 4) > pow(10, -4)) {pass = 0;} else {
  printf("Mg20 state #1 -> Mg20 state #1: 1-body J=0 T=1 sum rule: Pass\n");}

  fclose(in_file);

  in_file = fopen("../examples/output/mg20_mg20_1body_1_1.dens", "r");
  fscanf(in_file, "%d %d %d %d %d %d\n", &Jf, &Tf, &Ji, &Ti, &Jop, &Top);
  if (Jop != 0 || Top != 0) {printf("Error in density matrix\n"); exit(0);}
  mat = 0;
  for (int i = 0; i < 3; i++) {
   fscanf(in_file, "%d %d %d %d %f", &na, &ja, &nb, &jb, &dm);
   mat += sqrt(2.0)*sqrt(ja + 1.0)*dm;
  }
  mat *= 1.0/sqrt((Ji + 1.0)*(Ti + 1.0));
  if (fabs(mat - 4) > pow(10, -4)) {pass = 0;} else {
  printf("Mg20 state #2 -> Mg20 state #2: 1-body J=0 T=0 sum rule: Pass\n");}
  fscanf(in_file, "%*d\n");
  fscanf(in_file, "%d %d %d %d %d %d\n", &Jf, &Tf, &Ji, &Ti, &Jop, &Top);
  if (Jop != 0 || Top != 2) {printf("Error in density matrix\n"); exit(0);}
  mat = 0;
  for (int i = 0; i < 3; i++) {
   fscanf(in_file, "%d %d %d %d %f", &na, &ja, &nb, &jb, &dm);
   mat += sqrt(6.0)*sqrt(ja + 1.0)*dm;
  }
  mat *= 1.0/sqrt((Ji + 1.0)*(Ti + 1.0))*(sqrt(2.0/3.0));
  if (fabs(mat - 4) > pow(10, -4)) {pass = 0;} else {
  printf("Mg20 state #2 -> Mg20 state #2: 1-body J=0 T=1 sum rule: Pass\n");}

  fclose(in_file);

  in_file = fopen("../examples/output/mg20_mg20_1body_2_2.dens", "r");
  fscanf(in_file, "%d %d %d %d %d %d\n", &Jf, &Tf, &Ji, &Ti, &Jop, &Top);
  if (Jop != 0 || Top != 0) {printf("Error in density matrix\n"); exit(0);}
  mat = 0;
  for (int i = 0; i < 3; i++) {
   fscanf(in_file, "%d %d %d %d %f", &na, &ja, &nb, &jb, &dm);
   mat += sqrt(2.0)*sqrt(ja + 1.0)*dm;
  }
  mat *= 1.0/sqrt((Ji + 1.0)*(Ti + 1.0));
  if (fabs(mat - 4) > pow(10, -4)) {pass = 0;} else {
  printf("Mg20 state #3 -> Mg20 state #3: 1-body J=0 T=0 sum rule: Pass\n");}
  fscanf(in_file, "%*d\n");
  fscanf(in_file, "%d %d %d %d %d %d\n", &Jf, &Tf, &Ji, &Ti, &Jop, &Top);
  if (Jop != 0 || Top != 2) {printf("Error in density matrix\n"); exit(0);}
  mat = 0;
  for (int i = 0; i < 3; i++) {
   fscanf(in_file, "%d %d %d %d %f", &na, &ja, &nb, &jb, &dm);
   mat += sqrt(6.0)*sqrt(ja + 1.0)*dm;
  }
  mat *= 1.0/sqrt((Ji + 1.0)*(Ti + 1.0))*(sqrt(2.0/3.0));
  if (fabs(mat - 4) > pow(10, -4)) {pass = 0;} else {
  printf("Mg20 state #3 -> Mg20 state #3: 1-body J=0 T=1 sum rule: Pass\n");}
  fclose(in_file);
  printf("\n");

  sp = read_parameter_file("../examples/param_files/mg20_mg20_2body_test_T0.param");
  two_body_density(sp);
  free(sp);
  
  sp = read_parameter_file("../examples/param_files/mg20_mg20_2body_test_T1.param");
  two_body_density(sp);
  free(sp);
  
  mat = compute_2body_J0_T0_sum_rule("../examples/output/mg20_mg20_2body_J0_T0_0_0.dens")/sqrt(5.0);
  if (fabs(mat - 6) < pow(10, -4)) {printf("Mg20 state #1 -> Mg20 state #1: 2-body J=0 T=0 Sum Rule: Pass\n");}

  mat = compute_2body_J0_T1_sum_rule("../examples/output/mg20_mg20_2body_J0_T1_0_0.dens")*(sqrt(2.0/15.0));
  if (fabs(mat - 12) < pow(10, -4)) {printf("Mg20 state #1 -> Mg20 state #1: 2-body J=0 T=1 Sum Rule: Pass\n");}

  mat = compute_2body_J0_T0_sum_rule("../examples/output/mg20_mg20_2body_J0_T0_1_1.dens")/5.0;
  if (fabs(mat - 6) < pow(10, -4)) {printf("Mg20 state #2 -> Mg20 state #2: 2-body J=0 T=0 Sum Rule: Pass\n");}

  mat = compute_2body_J0_T1_sum_rule("../examples/output/mg20_mg20_2body_J0_T1_1_1.dens")*(sqrt(2.0/3.0)/5.0);
  if (fabs(mat - 12) < pow(10, -4)) {printf("Mg20 state #2 -> Mg20 state #2: 2-body J=0 T=1 Sum Rule: Pass\n");}

  mat = compute_2body_J0_T0_sum_rule("../examples/output/mg20_mg20_2body_J0_T0_2_2.dens")/(3.0*sqrt(5.0));
  if (fabs(mat - 6) < pow(10, -4)) {printf("Mg20 state #3 -> Mg20 state #3: 2-body J=0 T=0 Sum Rule: Pass\n");}

  mat = compute_2body_J0_T1_sum_rule("../examples/output/mg20_mg20_2body_J0_T1_2_2.dens")*(sqrt(2.0/15.0)/3.0);
  if (fabs(mat - 12) < pow(10, -4)) {printf("Mg20 state #3 -> Mg20 state #3: 2-body J=0 T=1 Sum Rule: Pass\n");}

  printf("\n");


  // Begin O20 -> O20 tests

  sp = read_parameter_file("../examples/param_files/o20_o20_1body_test.param");
  one_body_density(sp);

  in_file = fopen("../examples/output/o20_o20_1body_0_0.dens", "r");
  fscanf(in_file, "%d %d %d %d %d %d\n", &Jf, &Tf, &Ji, &Ti, &Jop, &Top);
  if (Jop != 0 || Top != 0) {printf("Error in density matrix\n"); exit(0);}
  mat = 0;
  for (int i = 0; i < 3; i++) {
   fscanf(in_file, "%d %d %d %d %f", &na, &ja, &nb, &jb, &dm);
   mat += sqrt(2.0)*sqrt(ja + 1.0)*dm;
  }
  mat *= 1.0/sqrt((Ji + 1.0)*(Ti + 1.0));
  if (fabs(mat - 4) > pow(10, -4)) {pass = 0;} else {
  printf("O20  state #1 -> O20  state #1: 1-body J=0 T=0 sum rule: Pass\n");}
  fscanf(in_file, "%*d\n");
  fscanf(in_file, "%d %d %d %d %d %d\n", &Jf, &Tf, &Ji, &Ti, &Jop, &Top);
  if (Jop != 0 || Top != 2) {printf("Error in density matrix\n"); exit(0);}
  mat = 0;
  for (int i = 0; i < 3; i++) {
   fscanf(in_file, "%d %d %d %d %f", &na, &ja, &nb, &jb, &dm);
   mat += sqrt(6.0)*sqrt(ja + 1.0)*dm;
  }
  mat *= 1.0/sqrt((Ji + 1.0)*(Ti + 1.0))*(-sqrt(2.0/3.0));
  if (fabs(mat + 4) > pow(10, -4)) {pass = 0;} else {
  printf("O20  state #1 -> O20  state #1: 1-body J=0 T=1 sum rule: Pass\n");}

  fclose(in_file);


  in_file = fopen("../examples/output/o20_o20_1body_1_1.dens", "r");
  fscanf(in_file, "%d %d %d %d %d %d\n", &Jf, &Tf, &Ji, &Ti, &Jop, &Top);
  if (Jop != 0 || Top != 0) {printf("Error in density matrix\n"); exit(0);}
  mat = 0;
  for (int i = 0; i < 3; i++) {
   fscanf(in_file, "%d %d %d %d %f", &na, &ja, &nb, &jb, &dm);
   mat += sqrt(2.0)*sqrt(ja + 1.0)*dm;
  }
  mat *= 1.0/sqrt((Ji + 1.0)*(Ti + 1.0));
  if (fabs(mat - 4) > pow(10, -4)) {pass = 0;} else {
  printf("O20  state #2 -> O20  state #2: 1-body J=0 T=0 sum rule: Pass\n");}
  fscanf(in_file, "%*d\n");
  fscanf(in_file, "%d %d %d %d %d %d\n", &Jf, &Tf, &Ji, &Ti, &Jop, &Top);
  if (Jop != 0 || Top != 2) {printf("Error in density matrix\n"); exit(0);}
  mat = 0;
  for (int i = 0; i < 3; i++) {
   fscanf(in_file, "%d %d %d %d %f", &na, &ja, &nb, &jb, &dm);
   mat += sqrt(6.0)*sqrt(ja + 1.0)*dm;
  }
  mat *= 1.0/sqrt((Ji + 1.0)*(Ti + 1.0))*(-sqrt(2.0/3.0));
  if (fabs(mat + 4) > pow(10, -4)) {pass = 0;} else {
  printf("O20  state #2 -> O20  state #2: 1-body J=0 T=1 sum rule: Pass\n");}

  fclose(in_file);

  in_file = fopen("../examples/output/o20_o20_1body_2_2.dens", "r");
  fscanf(in_file, "%d %d %d %d %d %d\n", &Jf, &Tf, &Ji, &Ti, &Jop, &Top);
  if (Jop != 0 || Top != 0) {printf("Error in density matrix\n"); exit(0);}
  mat = 0;
  for (int i = 0; i < 3; i++) {
   fscanf(in_file, "%d %d %d %d %f", &na, &ja, &nb, &jb, &dm);
   mat += sqrt(2.0)*sqrt(ja + 1.0)*dm;
  }
  mat *= 1.0/sqrt((Ji + 1.0)*(Ti + 1.0));
  if (fabs(mat - 4) > pow(10, -4)) {pass = 0;} else {
  printf("O20  state #3 -> O20  state #3: 1-body J=0 T=0 sum rule: Pass\n");}
  fscanf(in_file, "%*d\n");
  fscanf(in_file, "%d %d %d %d %d %d\n", &Jf, &Tf, &Ji, &Ti, &Jop, &Top);
  if (Jop != 0 || Top != 2) {printf("Error in density matrix\n"); exit(0);}
  mat = 0;
  for (int i = 0; i < 3; i++) {
   fscanf(in_file, "%d %d %d %d %f", &na, &ja, &nb, &jb, &dm);
   mat += sqrt(6.0)*sqrt(ja + 1.0)*dm;
  }
  mat *= 1.0/sqrt((Ji + 1.0)*(Ti + 1.0))*(-sqrt(2.0/3.0));
  if (fabs(mat + 4) > pow(10, -4)) {pass = 0;} else {
  printf("O20  state #3 -> O20  state #3: 1-body J=0 T=1 sum rule: Pass\n");}
  fclose(in_file);
  printf("\n");

  sp = read_parameter_file("../examples/param_files/o20_o20_2body_test_T0.param");
  two_body_density(sp);
  free(sp);
  
  sp = read_parameter_file("../examples/param_files/o20_o20_2body_test_T1.param");
  two_body_density(sp);
  free(sp);
  
  mat = compute_2body_J0_T0_sum_rule("../examples/output/o20_o20_2body_J0_T0_0_0.dens")/sqrt(5.0);
  if (fabs(mat - 6) < pow(10, -4)) {printf("O20  state #1 -> O20  state #1: 2-body J=0 T=0 Sum Rule: Pass\n");}

  mat = compute_2body_J0_T1_sum_rule("../examples/output/o20_o20_2body_J0_T1_0_0.dens")*(-sqrt(2.0/15.0));
  if (fabs(mat + 12) < pow(10, -4)) {printf("O20  state #1 -> O20  state #1: 2-body J=0 T=1 Sum Rule: Pass\n");}

  mat = compute_2body_J0_T0_sum_rule("../examples/output/o20_o20_2body_J0_T0_1_1.dens")/5.0;
  if (fabs(mat - 6) < pow(10, -4)) {printf("O20  state #2 -> O20  state #2: 2-body J=0 T=0 Sum Rule: Pass\n");}

  mat = compute_2body_J0_T1_sum_rule("../examples/output/o20_o20_2body_J0_T1_1_1.dens")*(-sqrt(2.0/3.0)/5.0);
  if (fabs(mat + 12) < pow(10, -4)) {printf("O20  state #2 -> O20  state #2: 2-body J=0 T=1 Sum Rule: Pass\n");}

  mat = compute_2body_J0_T0_sum_rule("../examples/output/o20_o20_2body_J0_T0_2_2.dens")/(3.0*sqrt(5.0));
  if (fabs(mat - 6) < pow(10, -4)) {printf("O20  state #3 -> O20  state #3: 2-body J=0 T=0 Sum Rule: Pass\n");}

  mat = compute_2body_J0_T1_sum_rule("../examples/output/o20_o20_2body_J0_T1_2_2.dens")*(-sqrt(2.0/15.0)/3.0);
  if (fabs(mat + 12) < pow(10, -4)) {printf("O20  state #3 -> O20  state #3: 2-body J=0 T=1 Sum Rule: Pass\n");}

  printf("\n");

  // Begin Na20 -> Ne20 Tests

  sp = read_parameter_file("../examples/param_files/na20_ne20_1body_test.param");
  one_body_density(sp);
  
  in_file = fopen("../examples/output/na20_ne20_1body_0_1.dens", "r");
  fscanf(in_file, "%d %d %d %d %d %d\n", &Jf, &Tf, &Ji, &Ti, &Jop, &Top);
  mat = 0;
  for (int i = 0; i < 3; i++) {
   fscanf(in_file, "%d %d %d %d %f", &na, &ja, &nb, &jb, &dm);
   mat += sqrt(2.0)*sqrt(ja + 1.0)*dm;
  }
  mat *= 1.0/sqrt((Ji + 1.0)*(Ti + 1.0));
  if (fabs(mat) > pow(10, -4)) {pass = 0;} else {
  printf("Na20 state #1 -> Ne20 state #2: 1-body J=0 T=1 sum rule: Pass\n");}
  in_file = fopen("../examples/output/na20_ne20_1body_2_2.dens", "r");
  fscanf(in_file, "%d %d %d %d %d %d\n", &Jf, &Tf, &Ji, &Ti, &Jop, &Top);
  mat = 0;
  for (int i = 0; i < 3; i++) {
   fscanf(in_file, "%d %d %d %d %f", &na, &ja, &nb, &jb, &dm);
   mat += sqrt(2.0)*sqrt(ja + 1.0)*dm;
  }
  mat *= 1.0/sqrt((Ji + 1.0)*(Ti + 1.0));
  if (fabs(mat) > pow(10, -4)) {pass = 0;} else {
  printf("Na20 state #3 -> Ne20 state #3: 1-body J=0 T=1 sum rule: Pass\n");}
  printf("\n");

  sp = read_parameter_file("../examples/param_files/na20_ne20_2body_test_T1.param");
  two_body_density(sp);

  mat = compute_2body_J0_T1_sum_rule("../examples/output/na20_ne20_2body_J0_T1_0_1.dens")/sqrt(5.0);

  if (fabs(mat) < pow(10, -4)) {printf("Na20 state #1 -> Ne20 state #2: 2-body J=0 T=1 Sum Rule: Pass\n");}

  mat = compute_2body_J0_T1_sum_rule("../examples/output/na20_ne20_2body_J0_T1_2_2.dens")/5.0;

  if (fabs(mat) < pow(10, -4)) {printf("Na20 state #3 -> Ne20 state #3: 2-body J=0 T=1 Sum Rule: Pass\n");}
  printf("\n");

  // Begin Ne20 -> Na20 Tests

  // One-body tests

  sp = read_parameter_file("../examples/param_files/ne20_na20_1body_test.param");
  one_body_density(sp);
  
  in_file = fopen("../examples/output/ne20_na20_1body_1_0.dens", "r");
  fscanf(in_file, "%d %d %d %d %d %d\n", &Jf, &Tf, &Ji, &Ti, &Jop, &Top);
  mat = 0;
  for (int i = 0; i < 3; i++) {
   fscanf(in_file, "%d %d %d %d %f", &na, &ja, &nb, &jb, &dm);
   mat += sqrt(2.0)*sqrt(ja + 1.0)*dm;
  }
  mat *= 1.0/sqrt((Ji + 1.0)*(Ti + 1.0));
  if (fabs(mat) > pow(10, -4)) {pass = 0;} else {
  printf("Ne20 state #2 -> Na20 state #1: 1-body J=0 T=1 sum rule: Pass\n");}
  in_file = fopen("../examples/output/ne20_na20_1body_2_2.dens", "r");
  fscanf(in_file, "%d %d %d %d %d %d\n", &Jf, &Tf, &Ji, &Ti, &Jop, &Top);
  mat = 0;
  for (int i = 0; i < 3; i++) {
   fscanf(in_file, "%d %d %d %d %f", &na, &ja, &nb, &jb, &dm);
   mat += sqrt(2.0)*sqrt(ja + 1.0)*dm;
  }
  mat *= 1.0/sqrt((Ji + 1.0)*(Ti + 1.0));
  if (fabs(mat) > pow(10, -4)) {pass = 0;} else {
  printf("Ne20 state #3 -> Na20 state #3: 1-body J=0 T=1 sum rule: Pass\n");}
  printf("\n");
  // Two-body tests

  sp = read_parameter_file("../examples/param_files/ne20_na20_2body_test_T1.param");
  two_body_density(sp);

  mat = compute_2body_J0_T1_sum_rule("../examples/output/ne20_na20_2body_J0_T1_1_0.dens")/sqrt(5.0);

  if (fabs(mat) < pow(10, -4)) {printf("Ne20 state #2 -> Na20 state #1: 2-body J=0 T=1 Sum Rule: Pass\n");}

  mat = compute_2body_J0_T1_sum_rule("../examples/output/ne20_na20_2body_J0_T1_2_2.dens")/5.0;

  if (fabs(mat) < pow(10, -4)) {printf("Ne20 state #3 -> Na20 state #3: 2-body J=0 T=1 Sum Rule: Pass\n");}
  printf("\n");


  // Begin F20 -> Ne20 Tests

  sp = read_parameter_file("../examples/param_files/f20_ne20_1body_test.param");
  one_body_density(sp);
  
  in_file = fopen("../examples/output/f20_ne20_1body_0_1.dens", "r");
  fscanf(in_file, "%d %d %d %d %d %d\n", &Jf, &Tf, &Ji, &Ti, &Jop, &Top);
  mat = 0;
  for (int i = 0; i < 3; i++) {
   fscanf(in_file, "%d %d %d %d %f", &na, &ja, &nb, &jb, &dm);
   mat += sqrt(2.0)*sqrt(ja + 1.0)*dm;
  }
  mat *= 1.0/sqrt((Ji + 1.0)*(Ti + 1.0));
  if (fabs(mat) > pow(10, -4)) {pass = 0;} else {
  printf("F20  state #1 -> Ne20 state #2: 1-body J=0 T=1 sum rule: Pass\n");}
  in_file = fopen("../examples/output/f20_ne20_1body_2_2.dens", "r");
  fscanf(in_file, "%d %d %d %d %d %d\n", &Jf, &Tf, &Ji, &Ti, &Jop, &Top);
  mat = 0;
  for (int i = 0; i < 3; i++) {
   fscanf(in_file, "%d %d %d %d %f", &na, &ja, &nb, &jb, &dm);
   mat += sqrt(2.0)*sqrt(ja + 1.0)*dm;
  }
  mat *= 1.0/sqrt((Ji + 1.0)*(Ti + 1.0));
  if (fabs(mat) > pow(10, -4)) {pass = 0;} else {
  printf("F20  state #3 -> Ne20 state #3: 1-body J=0 T=1 sum rule: Pass\n");}
  printf("\n");

  sp = read_parameter_file("../examples/param_files/f20_ne20_2body_test_T1.param");
  two_body_density(sp);

  mat = compute_2body_J0_T1_sum_rule("../examples/output/f20_ne20_2body_J0_T1_0_1.dens")/sqrt(5.0);

  if (fabs(mat) < pow(10, -4)) {printf("F20  state #1 -> Ne20 state #2: 2-body J=0 T=1 Sum Rule: Pass\n");}

  mat = compute_2body_J0_T1_sum_rule("../examples/output/f20_ne20_2body_J0_T1_2_2.dens")/5.0;

  if (fabs(mat) < pow(10, -4)) {printf("F20  state #3 -> Ne20 state #3: 2-body J=0 T=1 Sum Rule: Pass\n");}
  printf("\n");

  // Begin Ne20 -> F20 Tests

  // One-body tests

  sp = read_parameter_file("../examples/param_files/ne20_f20_1body_test.param");
  one_body_density(sp);
  
  in_file = fopen("../examples/output/ne20_f20_1body_1_0.dens", "r");
  fscanf(in_file, "%d %d %d %d %d %d\n", &Jf, &Tf, &Ji, &Ti, &Jop, &Top);
  mat = 0;
  for (int i = 0; i < 3; i++) {
   fscanf(in_file, "%d %d %d %d %f", &na, &ja, &nb, &jb, &dm);
   mat += sqrt(2.0)*sqrt(ja + 1.0)*dm;
  }
  mat *= 1.0/sqrt((Ji + 1.0)*(Ti + 1.0));
  if (fabs(mat) > pow(10, -4)) {pass = 0;} else {
  printf("Ne20 state #2 -> F20  state #1: 1-body J=0 T=1 sum rule: Pass\n");}
  in_file = fopen("../examples/output/ne20_f20_1body_2_2.dens", "r");
  fscanf(in_file, "%d %d %d %d %d %d\n", &Jf, &Tf, &Ji, &Ti, &Jop, &Top);
  mat = 0;
  for (int i = 0; i < 3; i++) {
   fscanf(in_file, "%d %d %d %d %f", &na, &ja, &nb, &jb, &dm);
   mat += sqrt(2.0)*sqrt(ja + 1.0)*dm;
  }
  mat *= 1.0/sqrt((Ji + 1.0)*(Ti + 1.0));
  if (fabs(mat) > pow(10, -4)) {pass = 0;} else {
  printf("Ne20 state #3 -> F20  state #3: 1-body J=0 T=1 sum rule: Pass\n");}
  printf("\n");
  // Two-body tests

  sp = read_parameter_file("../examples/param_files/ne20_f20_2body_test_T1.param");
  two_body_density(sp);

  mat = compute_2body_J0_T1_sum_rule("../examples/output/ne20_f20_2body_J0_T1_1_0.dens")/sqrt(5.0);

  if (fabs(mat) < pow(10, -4)) {printf("Ne20 state #2 -> F20  state #1: 2-body J=0 T=1 Sum Rule: Pass\n");}

  mat = compute_2body_J0_T1_sum_rule("../examples/output/ne20_f20_2body_J0_T1_2_2.dens")/5.0;

  if (fabs(mat) < pow(10, -4)) {printf("Ne20 state #3 -> F20  state #3: 2-body J=0 T=1 Sum Rule: Pass\n");}
  printf("\n");

  // Begin Mg20 -> Ne20 Tests
  
  sp = read_parameter_file("../examples/param_files/mg20_ne20_2body_test_T2.param");
  two_body_density(sp);

  mat = compute_2body_J0_T2_sum_rule("../examples/output/mg20_ne20_2body_J0_T2_0_0.dens")/sqrt(5.0);

  if (fabs(mat) < pow(10, -4)) {printf("Mg20 state #1 -> Ne20 state #1: 2-body J=0 T=2 Sum Rule: Pass\n");}

  mat = compute_2body_J0_T2_sum_rule("../examples/output/mg20_ne20_2body_J0_T2_1_1.dens")/5.0;

  if (fabs(mat) < pow(10, -4)) {printf("Mg20 state #2 -> Ne20 state #2: 2-body J=0 T=2 Sum Rule: Pass\n");}

  mat = compute_2body_J0_T2_sum_rule("../examples/output/mg20_ne20_2body_J0_T2_2_2.dens")/(3.0*sqrt(5.0));

  if (fabs(mat) < pow(10, -4)) {printf("Mg20 state #3 -> Ne20 state #3: 2-body J=0 T=0 Sum Rule: Pass\n");}
  printf("\n");
  free(sp);

  // Begin Ne20 -> Mg20 Tests
  
  sp = read_parameter_file("../examples/param_files/ne20_mg20_2body_test_T2.param");
  two_body_density(sp);

  mat = compute_2body_J0_T2_sum_rule("../examples/output/ne20_mg20_2body_J0_T2_0_0.dens")/sqrt(5.0);

  if (fabs(mat) < pow(10, -4)) {printf("Ne20 state #1 -> Mg20 state #1: 2-body J=0 T=2 Sum Rule: Pass\n");}

  mat = compute_2body_J0_T2_sum_rule("../examples/output/ne20_mg20_2body_J0_T2_1_1.dens")/5.0;

  if (fabs(mat) < pow(10, -4)) {printf("Ne20 state #2 -> Mg20 state #2: 2-body J=0 T=2 Sum Rule: Pass\n");}

  mat = compute_2body_J0_T2_sum_rule("../examples/output/ne20_mg20_2body_J0_T2_2_2.dens")/(3.0*sqrt(5.0));

  if (fabs(mat) < pow(10, -4)) {printf("Ne20 state #3 -> Mg20 state #3: 2-body J=0 T=0 Sum Rule: Pass\n");}
  printf("\n");
  free(sp);


  // Begin O20 -> Ne20 Tests
  
  sp = read_parameter_file("../examples/param_files/o20_ne20_2body_test_T2.param");
  two_body_density(sp);

  mat = compute_2body_J0_T2_sum_rule("../examples/output/o20_ne20_2body_J0_T2_0_0.dens")/sqrt(5.0);

  if (fabs(mat) < pow(10, -4)) {printf("O20  state #1 -> Ne20 state #1: 2-body J=0 T=2 Sum Rule: Pass\n");}

  mat = compute_2body_J0_T2_sum_rule("../examples/output/o20_ne20_2body_J0_T2_1_1.dens")/5.0;

  if (fabs(mat) < pow(10, -4)) {printf("O20  state #2 -> Ne20 state #2: 2-body J=0 T=2 Sum Rule: Pass\n");}

  mat = compute_2body_J0_T2_sum_rule("../examples/output/o20_ne20_2body_J0_T2_2_2.dens")/(3.0*sqrt(5.0));

  if (fabs(mat) < pow(10, -4)) {printf("O20  state #3 -> Ne20 state #3: 2-body J=0 T=0 Sum Rule: Pass\n");}
  printf("\n");
  free(sp);

  // Begin Ne20 -> O20 Tests
  
  sp = read_parameter_file("../examples/param_files/ne20_o20_2body_test_T2.param");
  two_body_density(sp);

  mat = compute_2body_J0_T2_sum_rule("../examples/output/ne20_o20_2body_J0_T2_0_0.dens")/sqrt(5.0);

  if (fabs(mat) < pow(10, -4)) {printf("Ne20 state #1 -> O20  state #1: 2-body J=0 T=2 Sum Rule: Pass\n");}

  mat = compute_2body_J0_T2_sum_rule("../examples/output/ne20_o20_2body_J0_T2_1_1.dens")/5.0;

  if (fabs(mat) < pow(10, -4)) {printf("Ne20 state #2 -> O20  state #2: 2-body J=0 T=2 Sum Rule: Pass\n");}

  mat = compute_2body_J0_T2_sum_rule("../examples/output/ne20_o20_2body_J0_T2_2_2.dens")/(3.0*sqrt(5.0));

  if (fabs(mat) < pow(10, -4)) {printf("Ne20 state #3 -> O20  state #3: 2-body J=0 T=0 Sum Rule: Pass\n");}
  printf("\n");
  free(sp);

  // Begin F20 -> Na20 Tests
  
  sp = read_parameter_file("../examples/param_files/f20_na20_2body_test_T2.param");
  two_body_density(sp);

  mat = compute_2body_J0_T2_sum_rule("../examples/output/f20_na20_2body_J0_T2_0_0.dens");
  if (fabs(mat - 5.0) < pow(10, -4)) {printf("F20  state #1 -> Na20 state #1: 2-body J=0 T=2 Sum Rule: Pass\n");}
  
  mat = compute_2body_J0_T2_sum_rule("../examples/output/f20_na20_2body_J0_T2_1_1.dens");
  if (fabs(mat - sqrt(35.0)) < pow(10, -4)) {printf("F20  state #2 -> Na20 state #2: 2-body J=0 T=2 Sum Rule: Pass\n");}

  mat = compute_2body_J0_T2_sum_rule("../examples/output/f20_na20_2body_J0_T2_2_2.dens");
  if (fabs(mat - 3.0*sqrt(5.0)) < pow(10, -4)) {printf("F20  state #3 -> Na20 state #3: 2-body J=0 T=0 Sum Rule: Pass\n");}

  printf("\n");
  free(sp);

  // Begin Na20 -> F20 Tests
  
  sp = read_parameter_file("../examples/param_files/na20_f20_2body_test_T2.param");
  two_body_density(sp);

  mat = compute_2body_J0_T2_sum_rule("../examples/output/na20_f20_2body_J0_T2_0_0.dens");
  if (fabs(mat - 5.0) < pow(10, -4)) {printf("Na20 state #1 -> F20  state #1: 2-body J=0 T=2 Sum Rule: Pass\n");}

  mat = compute_2body_J0_T2_sum_rule("../examples/output/na20_f20_2body_J0_T2_1_1.dens");
  if (fabs(mat - sqrt(35.0)) < pow(10, -4)) {printf("Na20 state #2 -> F20  state #2: 2-body J=0 T=2 Sum Rule: Pass\n");}

  mat = compute_2body_J0_T2_sum_rule("../examples/output/na20_f20_2body_J0_T2_2_2.dens");
  if (fabs(mat - 3.0*sqrt(5.0)) < pow(10, -4)) {printf("Na20 state #3 -> F20  state #3: 2-body J=0 T=0 Sum Rule: Pass\n");}

  printf("\n");
  free(sp);

  // Begin Al27 -> Si27 Tests

  // One-body tests

  sp = read_parameter_file("../examples/param_files/al27_si27_1body_test.param");
  one_body_density(sp);
  
  in_file = fopen("../examples/output/al27_si27_1body_0_0.dens", "r");
  fscanf(in_file, "%d %d %d %d %d %d\n", &Jf, &Tf, &Ji, &Ti, &Jop, &Top);
  mat = 0;
  for (int i = 0; i < 3; i++) {
   fscanf(in_file, "%d %d %d %d %f", &na, &ja, &nb, &jb, &dm);
   mat += sqrt(6.0)*sqrt(ja + 1.0)*dm;
  }
  mat *= 1.0/sqrt((Jf + 1.0)*(Tf + 1.0));
  if (fabs(mat + sqrt(3.0)) > pow(10, -4)) {pass = 0;} else {
  printf("Al27 state #1 -> Si27  state #1: 1-body J=0 T=1 sum rule: Pass\n");}

  in_file = fopen("../examples/output/al27_si27_1body_1_1.dens", "r");
  fscanf(in_file, "%d %d %d %d %d %d\n", &Jf, &Tf, &Ji, &Ti, &Jop, &Top);
  mat = 0;
  for (int i = 0; i < 3; i++) {
   fscanf(in_file, "%d %d %d %d %f", &na, &ja, &nb, &jb, &dm);
   mat += sqrt(6.0)*sqrt(ja + 1.0)*dm;
  }
  mat *= 1.0/sqrt((Ji + 1.0)*(Ti + 1.0));
  if (fabs(mat + sqrt(3.0)) > pow(10, -4)) {pass = 0;} else {
  printf("Al27 state #2 -> Si27  state #2: 1-body J=0 T=1 sum rule: Pass\n");}

  in_file = fopen("../examples/output/al27_si27_1body_2_2.dens", "r");
  fscanf(in_file, "%d %d %d %d %d %d\n", &Jf, &Tf, &Ji, &Ti, &Jop, &Top);
  mat = 0;
  for (int i = 0; i < 3; i++) {
   fscanf(in_file, "%d %d %d %d %f", &na, &ja, &nb, &jb, &dm);
   mat += sqrt(6.0)*sqrt(ja + 1.0)*dm;
  }
  mat *= 1.0/sqrt((Ji + 1.0)*(Ti + 1.0));
  if (fabs(mat + sqrt(3.0)) > pow(10, -4)) {pass = 0;} else {
  printf("Al27 state #3 -> Si27  state #3: 1-body J=0 T=1 sum rule: Pass\n");}

  printf("\n");
  free(sp);

  // Two-body tests
  sp = read_parameter_file("../examples/param_files/al27_si27_2body_test_T1.param");
  two_body_density(sp);

  mat = compute_2body_J0_T1_sum_rule("../examples/output/al27_si27_2body_J0_T1_0_0.dens")/6.0;
  if (fabs(mat - 10.0) < pow(10, -4)) {printf("Al27 state #1 -> Si27 state #1: 2-body J=0 T=2 Sum Rule: Pass\n");}

  mat = compute_2body_J0_T1_sum_rule("../examples/output/al27_si27_2body_J0_T1_1_1.dens")/(2*sqrt(3.0));
  if (fabs(mat - 10.0) < pow(10, -4)) {printf("Al27 state #2 -> Si27 state #2: 2-body J=0 T=2 Sum Rule: Pass\n");}

  mat = compute_2body_J0_T1_sum_rule("../examples/output/al27_si27_2body_J0_T1_2_2.dens")/(2*sqrt(6.0));
  if (fabs(mat - 10.0) < pow(10, -4)) {printf("Al27 state #3 -> Si27 state #3: 2-body J=0 T=0 Sum Rule: Pass\n");}
  
  printf("\n");
  free(sp);

  // Begin Si27 -> Al27 Tests

  // One-body tests

  sp = read_parameter_file("../examples/param_files/si27_al27_1body_test.param");
  one_body_density(sp);
  
  in_file = fopen("../examples/output/si27_al27_1body_0_0.dens", "r");
  fscanf(in_file, "%d %d %d %d %d %d\n", &Jf, &Tf, &Ji, &Ti, &Jop, &Top);
  mat = 0;
  for (int i = 0; i < 3; i++) {
   fscanf(in_file, "%d %d %d %d %f", &na, &ja, &nb, &jb, &dm);
   mat += sqrt(6.0)*sqrt(ja + 1.0)*dm;
  }
  mat *= 1.0/sqrt((Jf + 1.0)*(Tf + 1.0));
  if (fabs(mat + sqrt(3.0)) > pow(10, -4)) {pass = 0;} else {
  printf("Si27 state #1 -> Al27  state #1: 1-body J=0 T=1 sum rule: Pass\n");}

  in_file = fopen("../examples/output/si27_al27_1body_1_1.dens", "r");
  fscanf(in_file, "%d %d %d %d %d %d\n", &Jf, &Tf, &Ji, &Ti, &Jop, &Top);
  mat = 0;
  for (int i = 0; i < 3; i++) {
   fscanf(in_file, "%d %d %d %d %f", &na, &ja, &nb, &jb, &dm);
   mat += sqrt(6.0)*sqrt(ja + 1.0)*dm;
  }
  mat *= 1.0/sqrt((Ji + 1.0)*(Ti + 1.0));
  if (fabs(mat + sqrt(3.0)) > pow(10, -4)) {pass = 0;} else {
  printf("Si27 state #2 -> Al27  state #2: 1-body J=0 T=1 sum rule: Pass\n");}

  in_file = fopen("../examples/output/si27_al27_1body_2_2.dens", "r");
  fscanf(in_file, "%d %d %d %d %d %d\n", &Jf, &Tf, &Ji, &Ti, &Jop, &Top);
  mat = 0;
  for (int i = 0; i < 3; i++) {
   fscanf(in_file, "%d %d %d %d %f", &na, &ja, &nb, &jb, &dm);
   mat += sqrt(6.0)*sqrt(ja + 1.0)*dm;
  }
  mat *= 1.0/sqrt((Ji + 1.0)*(Ti + 1.0));
  if (fabs(mat + sqrt(3.0)) > pow(10, -4)) {pass = 0;} else {
  printf("Si27 state #3 -> Al27  state #3: 1-body J=0 T=1 sum rule: Pass\n");}

  printf("\n");
  free(sp);

  // Two-body tests
  sp = read_parameter_file("../examples/param_files/si27_al27_2body_test_T1.param");
  two_body_density(sp);

  mat = compute_2body_J0_T1_sum_rule("../examples/output/si27_al27_2body_J0_T1_0_0.dens")/6.0;
  if (fabs(mat - 10.0) < pow(10, -4)) {printf("Si27 state #1 -> Al27 state #1: 2-body J=0 T=2 Sum Rule: Pass\n");}

  mat = compute_2body_J0_T1_sum_rule("../examples/output/si27_al27_2body_J0_T1_1_1.dens")/(2*sqrt(3.0));
  if (fabs(mat - 10.0) < pow(10, -4)) {printf("Si27 state #2 -> Al27 state #2: 2-body J=0 T=2 Sum Rule: Pass\n");}

  mat = compute_2body_J0_T1_sum_rule("../examples/output/si27_al27_2body_J0_T1_2_2.dens")/(2*sqrt(6.0));
  if (fabs(mat - 10.0) < pow(10, -4)) {printf("Si27 state #3 -> Al27 state #3: 2-body J=0 T=0 Sum Rule: Pass\n");}
  
  printf("\n");
  free(sp);


  // Begin Al27 -> Mg27 Tests

  // One-body tests

  sp = read_parameter_file("../examples/param_files/al27_mg27_1body_test.param");
  one_body_density(sp);
  
  in_file = fopen("../examples/output/al27_mg27_1body_0_0.dens", "r");
  fscanf(in_file, "%d %d %d %d %d %d\n", &Jf, &Tf, &Ji, &Ti, &Jop, &Top);
  mat = 0;
  for (int i = 0; i < 3; i++) {
   fscanf(in_file, "%d %d %d %d %f", &na, &ja, &nb, &jb, &dm);
   mat += sqrt(6.0)*sqrt(ja + 1.0)*dm;
  }
  mat *= 1.0/sqrt((Jf + 1.0)*(Tf + 1.0));
  printf("Al Mat: %g\n", mat);
  if (fabs(mat) > pow(10, -4)) {pass = 0;} else {
  printf("Ne20 state #2 -> F20  state #1: 1-body J=0 T=1 sum rule: Pass\n");}
  in_file = fopen("../examples/output/al27_mg27_1body_0_1.dens", "r");
  fscanf(in_file, "%d %d %d %d %d %d\n", &Jf, &Tf, &Ji, &Ti, &Jop, &Top);
  mat = 0;
  for (int i = 0; i < 3; i++) {
   fscanf(in_file, "%d %d %d %d %f", &na, &ja, &nb, &jb, &dm);
   mat += sqrt(6.0)*sqrt(ja + 1.0)*dm;
  }
  mat *= 1.0/sqrt((Ji + 1.0)*(Ti + 1.0));
  printf("Al Mat: %g\n", mat);
  if (fabs(mat) > pow(10, -4)) {pass = 0;} else {
  printf("Ne20 state #3 -> F20  state #3: 1-body J=0 T=1 sum rule: Pass\n");}
  printf("\n");


  return pass;
}

double compute_2body_J0_T0_sum_rule(char *input_file) {
  FILE *in_file;
  in_file = fopen(input_file, "r");

  double mat = 0;
  float dm;
  int in1p, ij1p, in2p, ij2p, ij12p, it12p, in1, ij1, in2, ij2, ij12, it12;
  while(fscanf(in_file, "%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%f\n", &in1p, &ij1p, &in2p, &ij2p, &ij12p, &it12p, &in1, &ij1, &in2, &ij2, &ij12, &it12, &dm) == 13) {
    // The angular momentum are doubled in the file
    double j1 = ij1/2.0;
    double j2 = ij2/2.0;
    double j12 = ij12/2.0;
    double t12 = it12/2.0;
    double j1p = ij1p/2.0;
    double j2p = ij2p/2.0;
    double j12p = ij12p/2.0;
    double t12p = it12p/2.0;
    double m4 = 0.0;
    if (t12 != t12p) {continue;}
    if (j12 != j12p) {continue;}
    if ((in1 == in1p) && (j1 == j1p) && (in2 == in2p) && (j2 == j2p)) {
      m4 = 1.0;
    }
   
    if ((in1 == in2p) && (j1 == j2p) && (in2 == in1p) && (j2 == j1p)) {
      m4 += pow(-1.0, j1 + j2 + j12 + t12);
    }
    if (m4 == 0) {continue;}
    if ((in1 == in2) && (j1 == j2)) {m4 *= 1.0/sqrt(2.0);}
    if ((in1p == in2p) && (j1p == j2p)) {m4 *= 1.0/sqrt(2.0);}
    m4 *= sqrt(2.0*j12p + 1.0)*sqrt(2*t12 + 1);
    mat += m4*dm;
  }
  fclose(in_file);


  return mat;
}

double compute_2body_J0_T1_sum_rule(char *input_file) {
  FILE *in_file;
  in_file = fopen(input_file, "r");

  double mat = 0;
  float dm;
  int in1p, ij1p, in2p, ij2p, ij12p, it12p, in1, ij1, in2, ij2, ij12, it12;
  while(fscanf(in_file, "%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%f\n", &in1p, &ij1p, &in2p, &ij2p, &ij12p, &it12p, &in1, &ij1, &in2, &ij2, &ij12, &it12, &dm) == 13) {
    // The angular momentum are doubled in the file
    double j1 = ij1/2.0;
    double j2 = ij2/2.0;
    double j12 = ij12/2.0;
    double t12 = it12/2.0;
    double j1p = ij1p/2.0;
    double j2p = ij2p/2.0;
    double j12p = ij12p/2.0;
    double t12p = it12p/2.0;
    double m4 = 0.0;
    if (t12 != t12p) {continue;}
    if (j12 != j12p) {continue;}
    if (t12 != 1) {continue;}
    if ((in1 == in1p) && (j1 == j1p) && (in2 == in2p) && (j2 == j2p)) {
      m4 = 1.0;
    }
   
    if ((in1 == in2p) && (j1 == j2p) && (in2 == in1p) && (j2 == j1p)) {
      m4 += pow(-1.0, j1 + j2 + j12 + t12);
    }
    if (m4 == 0) {continue;}
    if ((in1 == in2) && (j1 == j2)) {m4 *= 1.0/sqrt(2.0);}
    if ((in1p == in2p) && (j1p == j2p)) {m4 *= 1.0/sqrt(2.0);}
    m4 *= sqrt(2.0*j12p + 1.0)*2.0*sqrt(6.0);
    mat += m4*dm;
  }
  fclose(in_file);

  

  return mat;
}

double compute_2body_J0_T2_sum_rule(char *input_file) {
  FILE *in_file;
  in_file = fopen(input_file, "r");

  double mat = 0;
  float dm;
  int in1p, ij1p, in2p, ij2p, ij12p, it12p, in1, ij1, in2, ij2, ij12, it12;
  while(fscanf(in_file, "%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%f\n", &in1p, &ij1p, &in2p, &ij2p, &ij12p, &it12p, &in1, &ij1, &in2, &ij2, &ij12, &it12, &dm) == 13) {
    // The angular momentum are doubled in the file
    double j1 = ij1/2.0;
    double j2 = ij2/2.0;
    double j12 = ij12/2.0;
    double t12 = it12/2.0;
    double j1p = ij1p/2.0;
    double j2p = ij2p/2.0;
    double j12p = ij12p/2.0;
    double t12p = it12p/2.0;
    double m4 = 0.0;
    if (t12 != t12p) {continue;}
    if (j12 != j12p) {continue;}
    if (t12 != 1) {continue;}
    if ((in1 == in1p) && (j1 == j1p) && (in2 == in2p) && (j2 == j2p)) {
      m4 = 1.0;
    }
   
    if ((in1 == in2p) && (j1 == j2p) && (in2 == in1p) && (j2 == j1p)) {
      m4 += pow(-1.0, j1 + j2 + j12 + t12);
    }
    if (m4 == 0) {continue;}
    if ((in1 == in2) && (j1 == j2)) {m4 *= 1.0/sqrt(2.0);}
    if ((in1p == in2p) && (j1p == j2p)) {m4 *= 1.0/sqrt(2.0);}
    m4 *= sqrt(2.0*j12p + 1.0)*sqrt(5.0);
    mat += m4*dm;
  }
  fclose(in_file);


  return mat;
}

