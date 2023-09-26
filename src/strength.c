#include "strength.h"

/* This file contains routines to generate one-body and two-body density matrices
   The input files are wavefunction files written by BIGSTICK
   The code relies on a factorization between proton and neutron Slater determinants
*/


void two_body_pivot(speedParams *sp) {
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
  wfnData *wd = read_binary_wfn_strength(wfn_file_initial, wfn_file_final, basis_file_initial, basis_file_final);

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

  if (wd->w_max_i > 0 || wd->w_max_f > 0) { 
    printf("Generating truncated jumps...\n");
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
      build_two_body_jumps_dmt0_trunc(wd->n_shells, wd->n_proton_i, mj_min_p_i, mj_max_p_i, num_mj_p_i, wd->n_sds_p_i, n_sds_p_int1, n_sds_p_int2, p1_array_f, p2_array_f, p0_list_i, p1_list_i, p2_list_i, wd->jz_shell, wd->l_shell, wd->w_shell, w_max_p_i);
      build_two_body_jumps_dmt0_trunc(wd->n_shells, wd->n_neutron_i, mj_min_n_i, mj_max_n_i, num_mj_n_i, wd->n_sds_n_i, n_sds_n_int1, n_sds_n_int2, n1_array_f, n2_array_f, n0_list_i, n1_list_i, n2_list_i, wd->jz_shell, wd->l_shell, wd->w_shell, w_max_n_i);
    } else if (mt_op == 2) {
      build_two_body_jumps_dmtp2_trunc(wd->n_shells, wd->n_proton_f, num_mj_p_i, mj_min_p_i, mj_max_p_i, num_mj_p_f, mj_min_p_f, mj_max_p_f, wd->n_sds_p_f, p2_list_f, wd->n_neutron_i, num_mj_n_i, mj_min_n_i, mj_max_n_i, num_mj_n_f, mj_min_n_f, mj_max_n_f, wd->n_sds_n_i, n2_list_i, wd->jz_shell, wd->l_shell, wd->w_shell, w_max_p_i, w_max_p_f, w_max_n_i, w_max_n_f);
    } else if (mt_op == -2) {
      build_two_body_jumps_dmtp2_trunc(wd->n_shells, wd->n_neutron_f, num_mj_n_i, mj_min_n_i, mj_max_n_i, num_mj_n_f, mj_min_n_f, mj_max_n_f, wd->n_sds_n_f, n2_list_f, wd->n_proton_i, num_mj_p_i, mj_min_p_i, mj_max_p_i, num_mj_p_f, mj_min_p_f, mj_max_p_f, wd->n_sds_p_i, p2_list_i, wd->jz_shell, wd->l_shell, wd->w_shell, w_max_n_i, w_max_n_f, w_max_p_i, w_max_p_f);
    } else if (mt_op == 1) {
      build_two_body_jumps_dmtp1_trunc(wd->n_shells, wd->n_proton_i, num_mj_p_i, mj_min_p_i, mj_max_p_i, wd->n_proton_f, num_mj_p_f, mj_min_p_f, mj_max_p_f, p1_list_i, p1_list_f, p2_array_f, wd->n_neutron_i, num_mj_n_i, mj_min_n_i, mj_max_n_i, wd->n_neutron_f, num_mj_n_f, mj_min_n_f, mj_max_n_f, n1_list_i, n2_list_i, n1_array_f, wd->jz_shell, wd->l_shell, wd->w_shell, w_max_p_i, w_max_p_f, w_max_n_i, w_max_n_f);    
    } else if (mt_op == -1) {
      build_two_body_jumps_dmtp1_trunc(wd->n_shells, wd->n_neutron_i, num_mj_n_i, mj_min_n_i, mj_max_n_i, wd->n_neutron_f, num_mj_n_f, mj_min_n_f, mj_max_n_f, n1_list_i, n1_list_f, n2_array_f, wd->n_proton_i, num_mj_p_i, mj_min_p_i, mj_max_p_i, wd->n_proton_f, num_mj_p_f, mj_min_p_f, mj_max_p_f, p1_list_i, p2_list_i, p1_array_f, wd->jz_shell, wd->l_shell, wd->w_shell, w_max_n_i, w_max_n_f, w_max_p_i, w_max_p_f);
    }
  } else {

  if (wd->same_basis) {
    // We assume that same basis implies that magnetic isospin is unchanged
      build_two_body_jumps_dmt0_alt(wd->n_shells, wd->n_proton_i, mj_min_p_i, mj_max_p_i, num_mj_p_i, wd->n_sds_p_i, n_sds_p_int1, n_sds_p_int2, p1_array_f, p2_array_f, p0_list_i, p1_list_i, p2_list_i, p2_list_if, wd->jz_shell, wd->l_shell);
      build_two_body_jumps_dmt0_alt(wd->n_shells, wd->n_neutron_i, mj_min_n_i, mj_max_n_i, num_mj_n_i, wd->n_sds_n_i, n_sds_n_int1, n_sds_n_int2, n1_array_f, n2_array_f, n0_list_i, n1_list_i, n2_list_i, n2_list_if, wd->jz_shell, wd->l_shell);
    } else if (mt_op == 2) {
      build_two_body_jumps_dmtp2(wd->n_shells, wd->n_proton_f, num_mj_p_i, mj_min_p_i, mj_max_p_i, num_mj_p_f, mj_min_p_f, mj_max_p_f, wd->n_sds_p_f, p2_list_f, wd->n_neutron_i, num_mj_n_i, mj_min_n_i, mj_max_n_i, num_mj_n_f, mj_min_n_f, mj_max_n_f, wd->n_sds_n_i, n2_list_i, wd->jz_shell, wd->l_shell);
    } else if (mt_op == -2) {
      build_two_body_jumps_dmtp2(wd->n_shells, wd->n_neutron_f, num_mj_n_i, mj_min_n_i, mj_max_n_i, num_mj_n_f, mj_min_n_f, mj_max_n_f, wd->n_sds_n_f, n2_list_f, wd->n_proton_i, num_mj_p_i, mj_min_p_i, mj_max_p_i, num_mj_p_f, mj_min_p_f, mj_max_p_f, wd->n_sds_p_i, p2_list_i, wd->jz_shell, wd->l_shell);
    } else if (mt_op == 1) {
      build_two_body_jumps_dmtp1(wd->n_shells, wd->n_proton_i, num_mj_p_i, mj_min_p_i, mj_max_p_i, wd->n_proton_f, num_mj_p_f, mj_min_p_f, mj_max_p_f, p1_list_i, p1_list_f, p2_array_f, wd->n_neutron_i, num_mj_n_i, mj_min_n_i, mj_max_n_i, wd->n_neutron_f, num_mj_n_f, mj_min_n_f, mj_max_n_f, n1_list_i, n2_list_i, n1_array_f, wd->jz_shell, wd->l_shell);    
    } else if (mt_op == -1) {
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
    strcpy(output_density_file, sp->out_file_base);
    sprintf(output_suffix, "_J%d_T%d_%d_%d.dens", j_op, t_op, psi_i, psi_f);
    strcat(output_density_file, output_suffix);
//    out_file = fopen(output_density_file, "w"); 
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
          printf("%d %d %d %d\n", i_orb1, i_orb2, i_orb3, i_orb4);
          if (pow(-1.0, wd->l_orb[i_orb1] + wd->l_orb[i_orb2] + wd->l_orb[i_orb3] + wd->l_orb[i_orb4]) != 1.0) {continue;} // Parity
          float j3 = wd->j_orb[i_orb4];
          // Allocate storage for each j12 and j34
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
                  if (mj1 + mj2 != mj3 + mj4) {continue;}
                  if (mt1 + mt2 - mt3 - mt4 != mt_op) {continue;}
                  for (int i = 0; i < sp->n_trans; i++) {density[i] = 0.0;}
                  if (mt_op == 2) {  // p^dag p^dag n n
                      trace_two_body_nodes_dmtp2_strength(a, b, d, c, num_mj_p_i, mj_min_p_i, num_mj_n_i, mj_min_n_i, n2_list_i, p2_list_f, wd, 0, sp->transition_list, density);

                  } else if (mt_op == -2) { // n^dag n^dag p p
                      trace_two_body_nodes_dmtp2_strength(a, b, d, c, num_mj_n_i, mj_min_n_i, num_mj_p_i, mj_min_p_i, p2_list_i, n2_list_f, wd, 1, sp->transition_list, density); 

		      //printf("%g\n", density[0]); 
                  } else {printf("Error in mt_op: %d\n", mt_op); exit(0);}
                  //printf("%g %g %g %g %g %g %g %g\n", mj1, mj2, mj3, mj4, mt1, mt2, mt3, mt4); 
                }
              }
            }
          }
        }
      }
    }
  } 
  int i_wfn = 0;
 
  for (unsigned int j = 0; j < wd->n_states_f; j++) {
              if(wd->bc_f[wd->n_eig_f*j] != 0.0) {//printf("%g\n", wd->bc_f[wd->n_eig_f*j]); 
		      i_wfn++;}
  }
  printf("num states: %d\n", i_wfn);

  modify_binary_wfn_strength("fe58_kbp.wfn", wd->bc_f); 

  return;
}


void trace_two_body_nodes_dmtp2_strength(int a, int b, int c, int d, int num_mj_p_i, float mj_min_p_i, int num_mj_n_i, float mj_min_n_i, sd_list** n2_list_i, sd_list** p2_list_f, wfnData* wd, int i_op, eigen_list* transition, double* density) {
/* 

*/
  int ns = wd->n_shells;
  double mat = double_gamow_teller2(a, b, c, d, wd->n_shell, wd->l_shell, wd->j_shell, wd->jz_shell);
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
            unsigned int index_i;
            unsigned int index_f;
	    int found_i = 0;
	    int found_f = 0;
            if (i_op == 0) {
              unsigned int p_hash_i = ppi + wd->n_sds_p_i*(pni % HASH_SIZE);
              wh_list* node3 = wd->wh_hash_i[p_hash_i];
              while (node3 != NULL) {
                if ((pni == node3->pn) && (ppi == node3->pp)) {
                  index_i = node3->index;
		  found_i = 1;
                  break;
                }
                node3 = node3->next;
              }
              if (found_i == 0) {node2 = node2->next; continue;}

              unsigned int p_hash_f = ppf + wd->n_sds_p_f*(pnf % HASH_SIZE);
              node3 = wd->wh_hash_f[p_hash_f];
              while (node3 != NULL) {
                if ((pnf == node3->pn) && (ppf == node3->pp)) {
                  index_f = node3->index;
		  found_f = 1;
                  break;
                }
                node3 = node3->next;
              }
              if (found_f == 0) {node2 = node2->next; continue;}
              eigen_list* eig_pair = transition;
              int i_trans = 0;
              while (eig_pair != NULL) {
                int psi_i = eig_pair->eig_i;
                int psi_f = eig_pair->eig_f;
	        wd->bc_f[psi_f + wd->n_eig_f*index_f] += mat*wd->bc_i[psi_i + wd->n_eig_i*index_i]*phase1*phase2;
                i_trans++;
                eig_pair = eig_pair->next;
              }
             } else {
              unsigned int p_hash_i = pni + wd->n_sds_p_i*(ppi % HASH_SIZE);
              wh_list* node3 = wd->wh_hash_i[p_hash_i];
              while (node3 != NULL) {
                if ((pni == node3->pp) && (ppi == node3->pn)) {
                  index_i = node3->index;
		  found_i = 1;
                  break;
                }
                node3 = node3->next;
              }
              if (found_i == 0) {node2 = node2->next; continue;}

              unsigned int p_hash_f = pnf + wd->n_sds_p_f*(ppf % HASH_SIZE);
              node3 = wd->wh_hash_f[p_hash_f];
              while (node3 != NULL) {
                if ((pnf == node3->pp) && (ppf == node3->pn)) {
                  index_f = node3->index;
		  found_f = 1;
                  break;
                }
                node3 = node3->next;
              }
              if (found_f == 0) {node2 = node2->next; continue;}
              eigen_list* eig_pair = transition;
              int i_trans = 0;
              while (eig_pair != NULL) {
                int psi_i = eig_pair->eig_i;
                int psi_f = eig_pair->eig_f;
	        wd->bc_f[psi_f + wd->n_eig_f*index_f] += mat*wd->bc_i[psi_i + wd->n_eig_i*index_i]*phase1*phase2;

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

double double_gamow_teller2(int a, int b, int c, int d, int* n_list, int* l_list, int* j_list, int* mj_list) {

//  if (a == c && b == d) {return 0.5;}
//  else {return 0;}

  int ija = j_list[a];
  int ijb = j_list[b];
  int ijc = j_list[c];
  int ijd = j_list[d];

  double ja = ija/2.0;
  double jb = ijb/2.0;
  double jc = ijc/2.0;
  double jd = ijd/2.0;

  int iNa = 2*(n_list[a] + 1) + l_list[a];
  int iNb = 2*(n_list[b] + 1) + l_list[b];
  int iNc = 2*(n_list[c] + 1) + l_list[c];
  int iNd = 2*(n_list[d] + 1) + l_list[d];
 
  double mja = mj_list[a]/2.0;
  double mjb = mj_list[b]/2.0;
  double mjc = mj_list[c]/2.0;
  double mjd = mj_list[d]/2.0;

  int ijab_max = (ija + ijb);
  int ijab_min = abs(ija - ijb);

  int ijcd_max = (ijc + ijd);
  int ijcd_min = abs(ijc - ijd);

  double mat = 0.0;

  for (int ijab = ijab_min; ijab <= ijab_max; ijab += 2) {
    double jab = ijab/2.0;
    if (fabs(mja + mjb) > jab) {continue;}
    for (int ijcd = ijcd_min; ijcd <= ijcd_max; ijcd += 2) {
      //if (ijab != ijcd) {continue;}
	double jcd = ijcd/2.0;
	if (fabs(mjc + mjd) > jcd) {continue;}
        double cg_fact = clebsch_gordan(ja, jb, jab, mja, mjb, mja + mjb)*clebsch_gordan(jc, jd, jcd, mjc, mjd, mjc + mjd);
	if (cg_fact == 0.0) {continue;}
        mat += cg_fact*compute_matrix_element_GT2(iNa, ija, iNb, ijb, ijab, iNc, ijc, iNd, ijd, ijcd, 2)/sqrt(2.0*jab + 1)*pow(-1.0, 2 - jcd + jab)*clebsch_gordan(2, jcd, jab, 0, mjc + mjd, mja + mjb);
    }
  }
  return mat;
}


double double_gamow_teller(int a, int b, int c, int d, int* n_list, int* l_list, int* j_list, int* mj_list) {

//  if (a == c && b == d) {return 0.5;}
//  else {return 0;}

  int ija = j_list[a];
  int ijb = j_list[b];
  int ijc = j_list[c];
  int ijd = j_list[d];

  double ja = ija/2.0;
  double jb = ijb/2.0;
  double jc = ijc/2.0;
  double jd = ijd/2.0;

  int iNa = 2*(n_list[a] + 1) + l_list[a];
  int iNb = 2*(n_list[b] + 1) + l_list[b];
  int iNc = 2*(n_list[c] + 1) + l_list[c];
  int iNd = 2*(n_list[d] + 1) + l_list[d];
 
  double mja = mj_list[a]/2.0;
  double mjb = mj_list[b]/2.0;
  double mjc = mj_list[c]/2.0;
  double mjd = mj_list[d]/2.0;

  int ijab_max = (ija + ijb);
  int ijab_min = abs(ija - ijb);

  int ijcd_max = (ijc + ijd);
  int ijcd_min = abs(ijc - ijd);

  double mat = 0.0;

  for (int ijab = ijab_min; ijab <= ijab_max; ijab += 2) {
    double jab = ijab/2.0;
    if (fabs(mja + mjb) > jab) {continue;}
    for (int ijcd = ijcd_min; ijcd <= ijcd_max; ijcd += 2) {
      if (ijab != ijcd) {continue;}
	double jcd = ijcd/2.0;
	if (fabs(mjc + mjd) > jcd) {continue;}
        double cg_fact = clebsch_gordan(ja, jb, jab, mja, mjb, mja + mjb)*clebsch_gordan(jc, jd, jcd, mjc, mjd, mjc + mjd);
	if (cg_fact == 0.0) {continue;}
        mat += cg_fact*compute_matrix_element_GT0(iNa, ija, iNb, ijb, ijab, iNc, ijc, iNd, ijd, ijcd, 2)/sqrt(2.0*jab + 1);
    }
  }
  return mat;
}


double double_fermi(int a, int b, int c, int d, int* n_list, int* l_list, int* j_list, int* mj_list) {

//  if (a == c && b == d) {return 0.5;}
//  else {return 0;}

  int ija = j_list[a];
  int ijb = j_list[b];
  int ijc = j_list[c];
  int ijd = j_list[d];

  double ja = ija/2.0;
  double jb = ijb/2.0;
  double jc = ijc/2.0;
  double jd = ijd/2.0;

  int iNa = 2*(n_list[a] + 1) + l_list[a];
  int iNb = 2*(n_list[b] + 1) + l_list[b];
  int iNc = 2*(n_list[c] + 1) + l_list[c];
  int iNd = 2*(n_list[d] + 1) + l_list[d];
 
  double mja = mj_list[a]/2.0;
  double mjb = mj_list[b]/2.0;
  double mjc = mj_list[c]/2.0;
  double mjd = mj_list[d]/2.0;

  int ijab_max = (ija + ijb);
  int ijab_min = abs(ija - ijb);

  int ijcd_max = (ijc + ijd);
  int ijcd_min = abs(ijc - ijd);

  double mat = 0.0;

  for (int ijab = ijab_min; ijab <= ijab_max; ijab += 2) {
    double jab = ijab/2.0;
    if (fabs(mja + mjb) > jab) {continue;}
    for (int ijcd = ijcd_min; ijcd <= ijcd_max; ijcd += 2) {
      if (ijab != ijcd) {continue;}
	double jcd = ijcd/2.0;
	if (fabs(mjc + mjd) > jcd) {continue;}
        double cg_fact = clebsch_gordan(ja, jb, jab, mja, mjb, mja + mjb)*clebsch_gordan(jc, jd, jcd, mjc, mjd, mjc + mjd);
	if (cg_fact == 0.0) {continue;}
        mat += cg_fact*compute_matrix_element_F0(iNa, ija, iNb, ijb, ijab, iNc, ijc, iNd, ijd, ijcd, 2)/sqrt(2.0*jab + 1);
    }
  }
  return mat;
}

double compute_matrix_element_GT2(int in1p, int ij1p, int in2p, int ij2p, int ij12p, int in1, int ij1, int in2, int ij2, int ij12, int it12) {

  double j1 = ij1/2.0;
  double j2 = ij2/2.0;
  double j12 = ij12/2.0;
  double t12 = it12/2.0;
  double j1p = ij1p/2.0;
  double j2p = ij2p/2.0;
  double j12p = ij12p/2.0;

  int l1, l2, l1p, l2p;
   
  l1p = get_l(in1p, ij1p);
  l2p = get_l(in2p, ij2p);
  l1 = get_l(in1, ij1);
  l2 = get_l(in2, ij2); 
    
  // The N's listed in the input file are energy quanta, we want radial quantum numbers
  double n1 = (in1 - l1)/2.0;
  double n2 = (in2 - l2)/2.0;
  double n1p = (in1p - l1p)/2.0;
  double n2p = (in2p - l2p)/2.0;
  double m4 = 0.0;
  // Lambda = Lambdap and S = SP
  for (int lambda = abs(l1 - l2); lambda <= (l1 + l2); lambda++) {
    if ((lambda < abs(l1p - l2p)) || (lambda > l1p + l2p)) {continue;} 
    int s_max = MIN(lambda + j12, 1);
    int s_min = abs(lambda - j12);
    for (int s = s_min; s <= s_max; s++) {
      // JJ -> LS coupling factors
      if (s != 1) {continue;}
      double fact = sqrt((2*lambda + 1)*(2*s + 1)*(2*j1 + 1)*(2*j2 + 1));
      fact *= sqrt((2*lambda + 1)*(2*s + 1)*(2*j1p + 1)*(2*j2p + 1));
      fact *= nine_j(l1, l2, lambda, 0.5, 0.5, s, j1, j2, j12);
      fact *= nine_j(l1p, l2p, lambda, 0.5, 0.5, s, j1p, j2p, j12p);
      // Un-reduce matrix element wrt total J
      fact *= sqrt(2*j12 + 1.0)*sqrt(2*j12p + 1.0)*pow(-1.0, lambda + 1 + j12p)*six_j(s, j12p, lambda, j12, s, 2);
      if (fact == 0.0) {continue;}
      // Un-reduced matrix element of sig(1) dot sig(2)
      double m1 = 2.0*sqrt(5.0);
      // Reduced matrix element of tau(1)_- tau(2)_-
//      m1 *= sqrt(5);
      double anti_symm = 0.0;
      if ((n1 == n1p) && (l1 == l1p) && (n2 == n2p) && (l2 == l2p)) {anti_symm = 1.0;}
      if ((n1 == n2p) && (l1 == l2p) && (n2 == n1p) && (l2 == l1p)) {anti_symm += pow(-1.0, t12 + l1 + l2 + lambda + s + 1);}
      if (anti_symm == 0) {continue;}
      m1 *= anti_symm;


      m1 *= fact;
      m4 += m1;
    }
  }

  if ((n1 == n2) && (j1 == j2) && (l1 == l2)) {m4 *= 1.0/2.0;}
  if ((n1p == n2p) && (j1p == j2p) && (l1p == l2p)) {m4 *= 1.0/2.0;}

//  printf("%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%g\n", in1p, ij1p, in2p, ij2p, ij12p, it12, in1, ij1, in2, ij2, ij12, it12, m4);

  return m4;
}


double compute_matrix_element_GT0(int in1p, int ij1p, int in2p, int ij2p, int ij12p, int in1, int ij1, int in2, int ij2, int ij12, int it12) {

  double j1 = ij1/2.0;
  double j2 = ij2/2.0;
  double j12 = ij12/2.0;
  double t12 = it12/2.0;
  double j1p = ij1p/2.0;
  double j2p = ij2p/2.0;
  double j12p = ij12p/2.0;

  int l1, l2, l1p, l2p;
   
  l1p = get_l(in1p, ij1p);
  l2p = get_l(in2p, ij2p);
  l1 = get_l(in1, ij1);
  l2 = get_l(in2, ij2); 
    
  // The N's listed in the input file are energy quanta, we want radial quantum numbers
  double n1 = (in1 - l1)/2.0;
  double n2 = (in2 - l2)/2.0;
  double n1p = (in1p - l1p)/2.0;
  double n2p = (in2p - l2p)/2.0;
  double m4 = 0.0;
  // Lambda = Lambdap and S = SP
  for (int lambda = abs(l1 - l2); lambda <= (l1 + l2); lambda++) {
    if ((lambda < abs(l1p - l2p)) || (lambda > l1p + l2p)) {continue;} 
    int s_max = MIN(lambda + j12, 1);
    int s_min = abs(lambda - j12);
    for (int s = s_min; s <= s_max; s++) {
      // JJ -> LS coupling factors
      double fact = sqrt((2*lambda + 1)*(2*s + 1)*(2*j1 + 1)*(2*j2 + 1));
      fact *= sqrt((2*lambda + 1)*(2*s + 1)*(2*j1p + 1)*(2*j2p + 1));
      fact *= nine_j(l1, l2, lambda, 0.5, 0.5, s, j1, j2, j12);
      fact *= nine_j(l1p, l2p, lambda, 0.5, 0.5, s, j1p, j2p, j12p);
      if (fact == 0.0) {continue;}
      // Un-reduce matrix element wrt total J
      fact *= sqrt(2*j12 + 1.0);
      // Un-reduced matrix element of sig(1) dot sig(2)
      double m1 = pow(-1.0, 1.0 + s)*six_j(s,0.5,0.5,1.0,0.5,0.5)*6.0;
      // Reduced matrix element of tau(1)_- tau(2)_-
//      m1 *= sqrt(5);
      double anti_symm = 0.0;
      if ((n1 == n1p) && (l1 == l1p) && (n2 == n2p) && (l2 == l2p)) {anti_symm = 1.0;}
      if ((n1 == n2p) && (l1 == l2p) && (n2 == n1p) && (l2 == l1p)) {anti_symm += pow(-1.0, t12 + l1 + l2 + lambda + s + 1);}
      if (anti_symm == 0) {continue;}
      m1 *= anti_symm;


      m1 *= fact;
      m4 += m1;
    }
  }

  if ((n1 == n2) && (j1 == j2) && (l1 == l2)) {m4 *= 1.0/2.0;}
  if ((n1p == n2p) && (j1p == j2p) && (l1p == l2p)) {m4 *= 1.0/2.0;}

//  printf("%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%g\n", in1p, ij1p, in2p, ij2p, ij12p, it12, in1, ij1, in2, ij2, ij12, it12, m4);

  return m4;
}


double compute_matrix_element_F0(int in1p, int ij1p, int in2p, int ij2p, int ij12p, int in1, int ij1, int in2, int ij2, int ij12, int it12) {

  double j1 = ij1/2.0;
  double j2 = ij2/2.0;
  double j12 = ij12/2.0;
  double t12 = it12/2.0;
  double j1p = ij1p/2.0;
  double j2p = ij2p/2.0;
  double j12p = ij12p/2.0;

  int l1, l2, l1p, l2p;
   
  l1p = get_l(in1p, ij1p);
  l2p = get_l(in2p, ij2p);
  l1 = get_l(in1, ij1);
  l2 = get_l(in2, ij2); 
    
  // The N's listed in the input file are energy quanta, we want radial quantum numbers
  double n1 = (in1 - l1)/2.0;
  double n2 = (in2 - l2)/2.0;
  double n1p = (in1p - l1p)/2.0;
  double n2p = (in2p - l2p)/2.0;
  double m4 = 0.0;
  // Lambda = Lambdap and S = SP
  for (int lambda = abs(l1 - l2); lambda <= (l1 + l2); lambda++) {
    if ((lambda < abs(l1p - l2p)) || (lambda > l1p + l2p)) {continue;} 
    int s_max = MIN(lambda + j12, 1);
    int s_min = abs(lambda - j12);
    for (int s = s_min; s <= s_max; s++) {
      // JJ -> LS coupling factors
      double fact = sqrt((2*lambda + 1)*(2*s + 1)*(2*j1 + 1)*(2*j2 + 1));
      fact *= sqrt((2*lambda + 1)*(2*s + 1)*(2*j1p + 1)*(2*j2p + 1));
      fact *= nine_j(l1, l2, lambda, 0.5, 0.5, s, j1, j2, j12);
      fact *= nine_j(l1p, l2p, lambda, 0.5, 0.5, s, j1p, j2p, j12p);
      if (fact == 0.0) {continue;}
      // Un-reduce matrix element wrt total J
      fact *= sqrt(2*j12 + 1.0);
      // Un-reduced matrix element of sig(1) dot sig(2)
      double m1 = 1.0;
      // Reduced matrix element of tau(1)_- tau(2)_-
//      m1 *= sqrt(5);
      double anti_symm = 0.0;
      if ((n1 == n1p) && (l1 == l1p) && (n2 == n2p) && (l2 == l2p)) {anti_symm = 1.0;}
//      if ((n1 == n2p) && (l1 == l2p) && (n2 == n1p) && (l2 == l1p)) {anti_symm += pow(-1.0, t12 + l1 + l2 + lambda + s + 1);}
      if (anti_symm == 0) {continue;}
      m1 *= anti_symm;


      m1 *= fact;
      m4 += m1;
    }
  }

//  if ((n1 == n2) && (j1 == j2) && (l1 == l2)) {m4 *= 1.0/sqrt(2.0);}
//  if ((n1p == n2p) && (j1p == j2p) && (l1p == l2p)) {m4 *= 1.0/sqrt(2.0);}

//  printf("%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%g\n", in1p, ij1p, in2p, ij2p, ij12p, it12, in1, ij1, in2, ij2, ij12, it12, m4);

  return m4;
}

int get_l(int N, int j) {
  int l;
  if ((N - (j + 1)/2) % 2 == 0) {
    l = (j + 1)/2;
  } else {
    l = (j - 1)/2;
  }
  return l;
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


void one_body_pivot(speedParams* sp) {
  // Reads in BIGSTICK basis/wavefunction (.trwfn) files along with
  // orbit definition (.sp) files and constructs the one-body density matrices
  // for each initial and final state eigenfunction
  // Initial and final wave functions must share the same orbit file
  // J_op and T_op are the total spin and isospin of the operator
  int iVerb = 1;
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
  wfnData *wd = read_binary_wfn_strength(wfn_file_initial, wfn_file_final, basis_file_initial, basis_file_final);

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
  
  if (wd->w_max_i > 0 || wd->w_max_f > 0) {if (iVerb) {printf("Truncation scheme detected...generating truncated jumps...\n");}
    int w_min_p_i = min_w(ns, wd->n_proton_i, wd->w_shell);
    int w_min_n_i = min_w(ns, wd->n_neutron_i, wd->w_shell);
    int w_min_p_f = min_w(ns, wd->n_proton_f, wd->w_shell);
    int w_min_n_f = min_w(ns, wd->n_neutron_f, wd->w_shell);
    
    int w_max_p_i = wd->w_max_i - w_min_n_i;
    int w_max_n_i = wd->w_max_i - w_min_p_i;
    int w_max_p_f = wd->w_max_f - w_min_n_f;
    int w_max_n_f = wd->w_max_f - w_min_p_f;


    if (wd->same_basis) {
      if (iVerb) {printf("Building initial and final state proton jumps...\n");}
      build_one_body_jumps_dmt0_trunc(wd->n_shells, wd->n_proton_i, mj_min_p_i, mj_max_p_i, num_mj_p_i, wd->n_sds_p_i, n_sds_p_int, p1_array_f, p0_list_i, p1_list_i, wd->jz_shell, wd->l_shell, wd->w_shell, w_max_p_i); 

      if (iVerb) {printf("Building initial and final state neutron jumps...\n");}
      build_one_body_jumps_dmt0_trunc(wd->n_shells, wd->n_neutron_i, mj_min_n_i, mj_max_n_i, num_mj_n_i, wd->n_sds_n_i, n_sds_n_int, n1_array_f, n0_list_i, n1_list_i, wd->jz_shell, wd->l_shell, wd->w_shell, w_max_n_i); 
    } else if (mt_op == 1) {
      if (iVerb) {printf("Delta MT = +1\n");}
      build_one_body_jumps_dmtp1_trunc(wd->n_shells, wd->n_proton_f, num_mj_p_i, mj_min_p_i, mj_max_p_i, num_mj_p_f, mj_min_p_f, mj_max_p_f, p1_list_f, wd->n_neutron_i, num_mj_n_i, mj_min_n_i, mj_max_n_i, n1_list_i, wd->jz_shell, wd->l_shell, wd->w_shell, w_max_p_i, w_max_p_f, w_max_n_i, w_max_n_f);
    } else if (mt_op == -1) {
      if (iVerb) {printf("Delta MT = -1\n");}
      build_one_body_jumps_dmtp1_trunc(wd->n_shells, wd->n_neutron_f, num_mj_n_i, mj_min_n_i, mj_max_n_i, num_mj_n_f, mj_min_n_f, mj_max_n_f, n1_list_f, wd->n_proton_i, num_mj_p_i, mj_min_p_i, mj_max_p_i, p1_list_i, wd->jz_shell, wd->l_shell, wd->w_shell, w_max_n_i, w_max_n_f, w_max_p_i, w_max_p_f);
    } else {printf("Error in delta MT\n"); exit(0);} 
  } else {
  if (wd->same_basis) {
      if (iVerb) {printf("Building initial and final state proton jumps...\n");}
      build_one_body_jumps_dmt0(wd->n_shells, wd->n_proton_i, mj_min_p_i, mj_max_p_i, num_mj_p_i, wd->n_sds_p_i, n_sds_p_int, p1_array_f, p0_list_i, p1_list_i, wd->jz_shell, wd->l_shell); 

      if (iVerb) {printf("Building initial and final state neutron jumps...\n");}
      build_one_body_jumps_dmt0(wd->n_shells, wd->n_neutron_i, mj_min_n_i, mj_max_n_i, num_mj_n_i, wd->n_sds_n_i, n_sds_n_int, n1_array_f, n0_list_i, n1_list_i, wd->jz_shell, wd->l_shell); 
    } else if (mt_op == 1) {
      if (iVerb) {printf("Delta MT = +1\n");}
      build_one_body_jumps_dmtp1(wd->n_shells, wd->n_proton_f, num_mj_p_i, mj_min_p_i, mj_max_p_i, num_mj_p_f, mj_min_p_f, mj_max_p_f, p1_list_f, wd->n_neutron_i, num_mj_n_i, mj_min_n_i, mj_max_n_i, n1_list_i, wd->jz_shell, wd->l_shell);
    } else if (mt_op == -1) {
      if (iVerb) {printf("Delta MT = -1\n");}
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
  //  printf("%d %d %g %g %g %g %g %g\n", psi_i, psi_f, ji, mji, jf, mjf, ti, tf);
    for (int j_op = (int) (fabs(jf - ji)); j_op <= (int) (ji + jf); j_op++) { 
      for (int t_op = MAX(0,abs(ti-tf)); t_op <= MIN(1,abs(ti+tf)); t_op++) {
	if (fabs(mt_op) > t_op) {continue;}
        cg_j = clebsch_gordan(j_op, ji, jf, 0, mji,mjf);
        if (fabs(cg_j) < pow(10, -6)) {if (iVerb) {printf("Warning: CG coefficient is zero for allowed transition: Ji = %g Jf = %g Jop = %d \nComputing the density matrix for this operator requires the raising/lowering operators (not currently supported.)\n This entry will be omitted even though there could be non-zero contributions.\n", ji, jf, j_op);} continue;}
        cg_t = clebsch_gordan(t_op, ti, tf, mt_op, mti, mtf);
        if (fabs(cg_t) < pow(10, -6)) {if (iVerb) {printf("Warning: CG coefficient is zero for allowed transition: Ti = %g Tf = %g Top = %d \nComputing the density matrix for this operator requires the raising/lowering operators (not currently supported.) \n This entry will be omitted even though there could be non-zero contributions.\n", ti, tf, t_op);} continue;}
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
	  printf("%g %g %g %g\n", j1, mj1, j2, mj2);
          if (mt1 - mt2 != mt_op) {continue;}

          for (int i = 0; i < sp->n_trans; i++) {
            density[i] = 0.0;
          }
          
          if ((mt1 == 0.5) && (mt2 == -0.5)) {
            trace_one_body_nodes_dmtp1_strength(a, b, num_mj_p_i, mj_min_p_i, num_mj_n_i, mj_min_n_i, n1_list_i, p1_list_f, wd, 0, sp->transition_list, density);
          } else {
            trace_one_body_nodes_dmtp1_strength(a, b, num_mj_n_i, mj_min_n_i, num_mj_p_i, mj_min_p_i, p1_list_i, n1_list_f, wd, 1, sp->transition_list, density);
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

  int i_wfn = 0;
 
  for (unsigned int j = 0; j < wd->n_states_f; j++) {
              if(wd->bc_f[wd->n_eig_f*j] != 0.0) {//printf("%g\n", wd->bc_f[wd->n_eig_f*j]); 
		      i_wfn++;}
  }
  printf("num states: %d\n", i_wfn);

  modify_binary_wfn_strength("al26_usdb.wfn", wd->bc_f); 


  return;
}

void trace_one_body_nodes_dmtp1_strength(int a, int b, int num_mj_p_i, float mj_min_p_i, int num_mj_n_i, float mj_min_n_i, sd_list** n1_list_i, sd_list** p1_list_f, wfnData* wd, int i_op, eigen_list *transition, double* density) {
/* 

*/
  double mat = gamow_teller(a, b, wd->n_shell, wd->l_shell, wd->j_shell, wd->jz_shell);

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
	      //density[i_trans] += wd->bc_i[psi_i + wd->n_eig_i*index_i]*wd->bc_f[psi_f + wd->n_eig_f*index_f]*phase1*phase2;
              unsigned int psi_index_i = psi_i + wd->n_eig_i*index_i;
              unsigned int psi_index_f = psi_f + wd->n_eig_f*index_f;
	      wd->bc_f[psi_f + wd->n_eig_f*index_f] += wd->bc_i[psi_index_i]*mat*phase1*phase2;

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
	      //density[i_trans] += wd->bc_i[psi_i + wd->n_eig_i*index_i]*wd->bc_f[psi_f + wd->n_eig_f*index_f]*phase1*phase2;
              unsigned int psi_index_i = psi_i + wd->n_eig_i*index_i;
	      unsigned int psi_index_f = psi_f + wd->n_eig_f*index_f;
	      wd->bc_f[psi_f + wd->n_eig_f*index_f] += wd->bc_i[psi_index_i]*mat*phase1*phase2;

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


double M_pn(int a, int b, int* n_list, int* l_list, int* j_list, int* mj_list) {


  int ija = j_list[a];
  int ijb = j_list[b];

  double ja = ija/2.0;
  double jb = ijb/2.0;

  int iNa = 2*(n_list[a] + 1) + l_list[a];
  int iNb = 2*(n_list[b] + 1) + l_list[b];

  int na = n_list[a];
  int nb = n_list[b];

  int la = l_list[a];
  int lb = l_list[b];

  double mja = mj_list[a]/2.0;
  double mjb = mj_list[b]/2.0;

  if (la != lb) {return 0.0;}
  if (iNa != iNb) {return 0.0;}
  if (fabs(mja - mjb) > 1) {return 0.0;}
  double cg_fact = pow(-1, 1 + ja - jb)*clebsch_gordan(1, jb, ja, mja - mjb, mjb, mja)/sqrt(2.0*ja + 1.0);  
  double mat = cg_fact*pow(-1, l)*6*sqrt(5.0)*sqrt((2.0*ja + 1.0)*(2.0*jb + 1.0)*(2.0*la + 1.0)*(2.0*lb + 1.0))*three_j(la, 2, lb, 0, 0, 0)*nine_j(la, lb, 2, 0.5, 0.5, 1, ja, jb, 1) ;
  
  return mat;
}


double gamow_teller(int a, int b, int* n_list, int* l_list, int* j_list, int* mj_list) {


  int ija = j_list[a];
  int ijb = j_list[b];

  double ja = ija/2.0;
  double jb = ijb/2.0;

  int iNa = 2*(n_list[a] + 1) + l_list[a];
  int iNb = 2*(n_list[b] + 1) + l_list[b];
 
  int la = l_list[a];
  int lb = l_list[b];

  double mja = mj_list[a]/2.0;
  double mjb = mj_list[b]/2.0;

  if (la != lb) {return 0.0;}
  if (iNa != iNb) {return 0.0;}
  if (fabs(mja - mjb) > 1) {return 0.0;}
  double cg_fact = pow(-1, 1 + ja - jb)*clebsch_gordan(1, jb, ja, mja - mjb, mjb, mja)/sqrt(2.0*ja + 1.0);  
  double mat = cg_fact*pow(-1, 1.5 + lb + ja)*sqrt(6.0)*sqrt((2.0*ja + 1.0)*(2.0*jb + 1.0))*six_j(0.5, ja, lb, jb, 0.5, 1);
  
  return mat;
}

