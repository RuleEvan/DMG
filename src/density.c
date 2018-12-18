#include "density.h"

/* This file contains routines to generate one-body and two-body density matrices
   The input files are wavefunction files written by BIGSTICK
   The code relies on a factorization between proton and neutron Slater determinants
*/

void two_body_density(int j_op, int t_op) {

  // Read in data 
  char wfn_file_initial[] = "ge64_basis.wfn";
  char wfn_file_final[] = "se64_basis.wfn";
  char basis_file_initial[] = "ge64_basis.bas";
  char basis_file_final[] = "se64_basis.bas";
  char orbit_file[] = "GCN2850.sps";
  wfnData *wd = read_binary_wfn_data(wfn_file_initial, wfn_file_final, basis_file_initial, basis_file_final, orbit_file);
//  wfnData *wd = read_wfn_data(wfn_file_initial, wfn_file_final, orbit_file);
  float* j_store = (float*) malloc(sizeof(float)*4);

  int ns = wd->n_shells;

  // Determine number of intermediate Slater determinants
  // after the action of one or two annihilation operators
  int n_sds_p_int1 = get_num_sds(wd->n_shells, wd->n_proton_f - 1);
  int n_sds_p_int2 = get_num_sds(wd->n_shells, wd->n_proton_f - 2);
  int n_sds_n_int1 = get_num_sds(wd->n_shells, wd->n_neutron_f - 1);
  int n_sds_n_int2 = get_num_sds(wd->n_shells, wd->n_neutron_f - 2);

  printf("Intermediate SDs: p1: %d n1: %d, p2: %d, n2: %d\n", n_sds_p_int1, n_sds_n_int1, n_sds_p_int2, n_sds_n_int2);

  // Determine the min and max mj values of proton/neutron SDs
  float mj_min_p_i = min_mj(ns, wd->n_proton_i, wd->jz_shell);
  float mj_max_p_i = max_mj(ns, wd->n_proton_i, wd->jz_shell);
  float mj_min_n_i = min_mj(ns, wd->n_neutron_i, wd->jz_shell);
  float mj_max_n_i = max_mj(ns, wd->n_neutron_i, wd->jz_shell);

  // Account for the factor that m_p + m_n = 0 in determining min/max mj
  mj_min_p_i = MAX(mj_min_p_i, -mj_max_n_i);
  mj_max_p_i = MIN(mj_max_p_i, -mj_min_n_i);
  mj_min_n_i = MAX(mj_min_n_i, -mj_max_p_i);
  mj_max_n_i = MIN(mj_max_n_i, -mj_min_p_i);

  // Determine number of sectors (mp, mn)
  int num_mj_i = mj_max_p_i - mj_min_p_i + 1;
  printf("Number of m_j sectors:%d\n", num_mj_i);

  // Allocate space for jump lists
  wf_list **p0_list_i = (wf_list**) calloc(2*num_mj_i, sizeof(wf_list*));
  wf_list **n0_list_i = (wf_list**) calloc(2*num_mj_i, sizeof(wf_list*));
  sd_list **p1_list_i = (sd_list**) calloc(2*ns*num_mj_i, sizeof(sd_list*));
  sd_list **n1_list_i = (sd_list**) calloc(2*ns*num_mj_i, sizeof(sd_list*));
  sd_list **p2_list_i = (sd_list**) calloc(2*ns*ns*num_mj_i, sizeof(sd_list*));
  sd_list **n2_list_i = (sd_list**) calloc(2*ns*ns*num_mj_i, sizeof(sd_list*));
  sd_list **p2_list_f = (sd_list**) calloc(2*ns*ns*num_mj_i, sizeof(sd_list*));
  sd_list **n2_list_f = (sd_list**) calloc(2*ns*ns*num_mj_i, sizeof(sd_list*));
  int* p1_array_f = (int*) calloc(ns*n_sds_p_int1, sizeof(int));
  int* p2_array_f = (int*) calloc(ns*ns*n_sds_p_int2, sizeof(int));
  int* n1_array_f = (int*) calloc(ns*n_sds_n_int1, sizeof(int));
  int* n2_array_f = (int*) calloc(ns*ns*n_sds_n_int2, sizeof(int));

  // Determine one/two-body jumps, using special routine if initial and final bases are the same
  if (wd->same_basis) {
    printf("Building proton jumps...\n");
    build_two_body_jumps_i_and_f(wd->n_shells, wd->n_proton_i, mj_min_p_i, mj_max_p_i, num_mj_i, wd->n_sds_p_i, n_sds_p_int1, n_sds_p_int2, p1_array_f, p2_array_f, p0_list_i, p1_list_i, p2_list_i, wd->jz_shell, wd->l_shell);
    printf("Done.\n");
    printf("Building neutron jumps...\n");
    build_two_body_jumps_i_and_f(wd->n_shells, wd->n_neutron_i, mj_min_n_i, mj_max_n_i, num_mj_i, wd->n_sds_n_i, n_sds_n_int1, n_sds_n_int2, n1_array_f, n2_array_f, n0_list_i, n1_list_i, n2_list_i, wd->jz_shell, wd->l_shell);
    printf("Done.\n");
  } else {
    printf("Building initial state proton jumps...\n");
    build_two_body_jumps_i(wd->n_shells, wd->n_proton_i, mj_min_p_i, mj_max_p_i, num_mj_i, wd->n_sds_p_i, p0_list_i, p1_list_i, p2_list_i, wd->jz_shell, wd->l_shell);
    printf("Done\n");
    printf("Building final state proton jumps...\n");
    build_two_body_jumps_f(wd->n_shells, wd->n_proton_f, mj_min_p_i, mj_max_p_i, num_mj_i, wd->n_sds_p_f, n_sds_p_int1, n_sds_p_int2, p1_array_f, p2_array_f, p2_list_f, wd->jz_shell, wd->l_shell);
    printf("Done\n");
    printf("Building initial state neutron jumps...\n");
    build_two_body_jumps_i(wd->n_shells, wd->n_neutron_i, mj_min_n_i, mj_max_n_i, num_mj_i, wd->n_sds_n_i, n0_list_i, n1_list_i, n2_list_i, wd->jz_shell, wd->l_shell);
    printf("Building final state neutron jumps...\n");
    build_two_body_jumps_f(wd->n_shells, wd->n_neutron_f, mj_min_n_i, mj_max_n_i, num_mj_i, wd->n_sds_n_f, n_sds_n_int1, n_sds_n_int2, n1_array_f, n2_array_f, n2_list_f, wd->jz_shell, wd->l_shell);
    printf("Done.\n");
  }

  FILE *out_file;
  out_file = fopen("ne-mg_fermi_density", "w"); 
  for (int psi_i = 0; psi_i < 1; psi_i++) {
    float ji = wd->j_nuc_i[psi_i];
    float ti = wd->t_nuc_i[psi_i];
    // Many-body states do not have mt = 0
    float mti = 0.5*(wd->n_proton_i - wd->n_neutron_i);
    // Loop over final many-body wave functions
    for (int psi_f = 0; psi_f < 1; psi_f++) {
      float cg_j = 0.0;
      float cg_t = 0.0;
      float jf = wd->j_nuc_f[psi_f];
      float tf = wd->t_nuc_f[psi_f];
      // Many-body states do not have mt = 0
      float mtf = 0.5*(wd->n_proton_f - wd->n_neutron_f);
      float mt_op = mtf - mti;
      // CG factors for coupling initial and final states to operator
      cg_j = clebsch_gordan(j_op, ji, jf, 0, 0, 0);
      if (cg_j == 0.0) {continue;}
      cg_t = clebsch_gordan(t_op, ti, tf, mt_op, mti, mtf);
      if (cg_t == 0.0) {continue;}
      cg_j *= pow(-1.0, jf - ji)/sqrt(2*jf + 1);
      cg_t *= pow(-1.0, tf - ti)/sqrt(2*tf + 1);
      printf("Initial state: #%d J: %g T: %g MT: %g Final State: #%d J: %g T: %g MT: %g\n", psi_i + 1, ji, ti, mti, psi_f + 1, jf, tf, mtf);
      float mat_test = 0.0;
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
              j_store = realloc(j_store, sizeof(float)*4*j_dim);
              for (int k = 0; k < 4*j_dim; k++) {
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
                      float d1 = 0.0;
                      if (pow(-1.0, wd->l_shell[a] + wd->l_shell[b] + wd->l_shell[c] + wd->l_shell[d]) != 1.0) {continue;} // Parity
                      if ((mt3 == 0.5) && (mt4 == 0.5)) { // 
                        if ((mt1 == 0.5) && (mt2 == 0.5)) { // 2 proton creation operators + 2 proton annihilation operators
                          d1 = trace_a4_nodes(a, b, c, d, num_mj_i, n_sds_p_int2, p2_array_f, p2_list_i, n0_list_i, wd, psi_i, psi_f, 0);
                        } else if ((mt1 == -0.5) && (mt2 == -0.5)) { // 2 neutron creation operators and two proton ann. operators
                          if ((fabs(mj1 + mj2) > (mj_max_n_i - mj_min_n_i)) || (fabs(mj3 + mj4) > (mj_max_p_i - mj_min_p_i))) {printf("Saved time\n"); continue;}
                          d1 = trace_a22_nodes(a, b, c, d, num_mj_i, p2_list_i, n2_list_f, wd, psi_i, psi_f, 0);
                        }  
                      } else if ((mt3 == -0.5) && (mt4 == -0.5)) {
                        if ((mt1 == -0.5) && (mt2 == -0.5)) { //2 n cr. and 2 n ann. operators
                          d1 = trace_a4_nodes(a, b, c, d, num_mj_i, n_sds_n_int2, n2_array_f, n2_list_i, p0_list_i, wd, psi_i, psi_f, 1);
                        } else if ((mt1 == 0.5) && (mt2 == 0.5)) {// 2 p cr. and 2 n ann. operators
                          if ((fabs(mj1 + mj2) > (mj_max_p_i - mj_min_p_i)) || (fabs(mj3 + mj4) > (mj_max_n_i - mj_min_n_i))) {printf("Saved time\n");continue;}
                          d1 = trace_a22_nodes(a, b, c, d, num_mj_i, n2_list_i, p2_list_f, wd, psi_i, psi_f, 1);
                        }
                      } else if ((mt1 == 0.5) && (mt2 == -0.5) && (mt3 == 0.5) && (mt4 == -0.5)) {
                          d1 = -trace_a20_nodes(a, d, b, c, num_mj_i, n_sds_p_int1, n_sds_n_int1, p1_list_i, n1_list_i, p1_array_f, n1_array_f, wd, psi_i, psi_f);
                      } else if ((mt1 == -0.5) && (mt2 == 0.5) && (mt3 == 0.5) && (mt4 == -0.5)) {
                          d1 = trace_a20_nodes(b, d, a, c, num_mj_i, n_sds_p_int1, n_sds_n_int1, p1_list_i, n1_list_i, p1_array_f, n1_array_f, wd, psi_i, psi_f);
                      } else if ((mt1 == 0.5) && (mt2 == -0.5) && (mt3 == -0.5) && (mt4 == 0.5)) {
                          d1 = trace_a20_nodes(a, c, b, d, num_mj_i, n_sds_p_int1, n_sds_n_int1,p1_list_i, n1_list_i, p1_array_f, n1_array_f, wd, psi_i, psi_f);
                      } else if ((mt1 == -0.5) && (mt2 == 0.5) && (mt3 == -0.5) && (mt4 == 0.5)) {
                          d1 = -trace_a20_nodes(b, c, a, d, num_mj_i, n_sds_p_int1, n_sds_n_int1, p1_list_i, n1_list_i, p1_array_f, n1_array_f, wd, psi_i, psi_f);
                      }
                      printf("%d %d %d %d %g\n", a, b, c, d, d1);
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
                              float d2 = cg_j12*cg_j34*cg_jop;
                              d2 *= cg_t12*cg_t34*cg_top;
                              d2 *= pow(-1.0, -j34 + j12 - t34 + t12)/sqrt((2*j12 + 1)*(2*t12 + 1));
                              if (d2 == 0.0) {continue;}
                              j_store[t12 + 2*(t34 + 2*((j12 - j_min_12)*j_dim_34 + (j34 - j_min_34)))] += d1*d2/(cg_j*cg_t);
                              if ((i_orb1 == i_orb2) && (mt1 == mt2)) {
                                j_store[t12 + 2*(t34 + 2*((j12 - j_min_12)*j_dim_34 + (j34 - j_min_34)))] += pow(-1.0, j1 + j2 - j12 - t12)*d1*d2/(cg_j*cg_t);
                              }
                              if ((i_orb3 == i_orb4) && (mt3 == mt4)) {
                                j_store[t12 + 2*(t34 + 2*((j12 - j_min_12)*j_dim_34 + (j34 - j_min_34)))] += pow(-1.0, j3 + j4 - j34 - t34)*d1*d2/(cg_j*cg_t);
                              }
                              if ((i_orb1 == i_orb2) && (mt1 == mt2) && (i_orb3 == i_orb4) && (mt3 == mt4)) {
                                j_store[t12 + 2*(t34 + 2*((j12 - j_min_12)*j_dim_34 + (j34 - j_min_34)))] += pow(-1.0, j1 + j2 - j12 - t12 + j3 + j4 - j34 - t34)*d1*d2/(cg_j*cg_t);
                              }

                            }
                          }
                        }
                      }
                    }
                  }
                }
              }
           
              for (int t12 = 0; t12 <= 1; t12++) {
                for (int t34 = 0; t34 <= 1; t34++) {
                  for (int ij12 = 0; ij12 < j_dim_12; ij12++) {
                    for (int ij34 = 0; ij34 < j_dim_34; ij34++) {
                      if (fabs(j_store[t12 + 2*(t34 + 2*(ij12*j_dim_34 + ij34))]) < pow(10, -8)) {continue;}
                      fprintf(out_file, "%d,%g,%d,%g,%d,%d,%d,%g,%d,%g,%d,%d,%g\n", 2*wd->n_orb[i_orb1], 2*wd->j_orb[i_orb1], 2*wd->n_orb[i_orb2], 2*wd->j_orb[i_orb2], 2*(ij12 + j_min_12), 2*t12, 2*wd->n_orb[i_orb3], 2*wd->j_orb[i_orb3], 2*wd->n_orb[i_orb4], 2*wd->j_orb[i_orb4], 2*(ij34 + j_min_34), 2*t34, j_store[t12 + 2*(t34 + 2*(ij12*j_dim_34 + ij34))]);
                      if (i_orb1 != i_orb2) { 
                        fprintf(out_file, "%d,%g,%d,%g,%d,%d,%d,%g,%d,%g,%d,%d,%g\n", 2*wd->n_orb[i_orb2], 2*wd->j_orb[i_orb2], 2*wd->n_orb[i_orb1], 2*wd->j_orb[i_orb1], 2*(ij12 + j_min_12), 2*t12, 2*wd->n_orb[i_orb3], 2*wd->j_orb[i_orb3], 2*wd->n_orb[i_orb4], 2*wd->j_orb[i_orb4], 2*(ij34 + j_min_34), 2*t34, pow(-1.0, j1 + j2 - ij12 - j_min_12 - t12)*j_store[t12 + 2*(t34 + 2*(ij12*j_dim_34 + ij34))]);
                      }
                      if (i_orb3 != i_orb4) {
                        fprintf(out_file, "%d,%g,%d,%g,%d,%d,%d,%g,%d,%g,%d,%d,%g\n", 2*wd->n_orb[i_orb1], 2*wd->j_orb[i_orb1], 2*wd->n_orb[i_orb2], 2*wd->j_orb[i_orb2], 2*(ij12 + j_min_12), 2*t12, 2*wd->n_orb[i_orb4], 2*wd->j_orb[i_orb4], 2*wd->n_orb[i_orb3], 2*wd->j_orb[i_orb3], 2*(ij34 + j_min_34), 2*t34, pow(-1.0, j3 + j4 - ij34 - j_min_34 - t34)*j_store[t12 + 2*(t34 + 2*(ij12*j_dim_34 + ij34))]);
                      }
                      if ((i_orb1 != i_orb2) && (i_orb3 != i_orb4)) {
                        fprintf(out_file, "%d,%g,%d,%g,%d,%d,%d,%g,%d,%g,%d,%d,%g\n", 2*wd->n_orb[i_orb2], 2*wd->j_orb[i_orb2], 2*wd->n_orb[i_orb1], 2*wd->j_orb[i_orb1], 2*(ij12 + j_min_12), 2*t12, 2*wd->n_orb[i_orb4], 2*wd->j_orb[i_orb4], 2*wd->n_orb[i_orb3], 2*wd->j_orb[i_orb3], 2*(ij34 + j_min_34), 2*t34, pow(-1.0, j1 + j2 + j3 + j4 - ij12 - j_min_12 - ij34 - j_min_34 - t12 - t34)*j_store[t12 + 2*(t34 + 2*(ij12*j_dim_34 + ij34))]);

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
  } 
  free(j_store); 
  fclose(out_file);
  return;
}

float trace_a4_nodes(int a, int b, int c, int d, int num_mj, int n_sds_int2, int* p2_array_f, sd_list** p2_list_i, wf_list** n0_list_i, wfnData* wd, int psi_i, int psi_f, int i_op) {
  float total = 0.0;
  int ns = wd->n_shells;
  for (int imj = 0; imj < num_mj; imj++) {
    sd_list* node1 = p2_list_i[imj + num_mj*(c + d*ns)];
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
      wf_list* node2 = n0_list_i[num_mj - imj - 1];
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
          total += wd->bc_i[psi_i + wd->n_eig_i*index_i]*wd->bc_f[psi_f + wd->n_eig_f*index_f]*phase1*phase2;
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
          total += wd->bc_i[psi_i + wd->n_eig_i*index_i]*wd->bc_f[psi_f + wd->n_eig_f*index_f]*phase1*phase2;

        } 
        node2 = node2->next;
      }
    }
  }
  return total;
}

float trace_a22_nodes(int a, int b, int c, int d, int num_mj, sd_list** a2_list_i, sd_list** a2_list_f, wfnData* wd, int psi_i, int psi_f, int i_op) {
/* 

*/
  float total = 0.0;
  int ns = wd->n_shells;
  for (int ipar = 0; ipar <=1; ipar++) {
    for (int imj = 0; imj < num_mj; imj++) {
      sd_list* node1 = a2_list_i[ipar + 2*(imj + num_mj*(c + d*ns))];
      // Loop over final states resulting from 2x a_op
      while (node1 != NULL) {
        unsigned int ppf = node1->pn; // Get final state p_f
        unsigned int ppi = node1->pi; // Get initial state p_i
        int phase1 = node1->phase;
        // Get list of n_i associated to p_i
        sd_list* node2 = a2_list_f[ipar + 2*(num_mj - imj - 1 + num_mj*(a + b*ns))]; //hash corresponds to a_op operators
        // Loop over n_f  
        //int m_pf = m_from_p(ppf, ns, npp, wd->jz_shell);
        while (node2 != NULL) {
          unsigned int pni = node2->pi;
          unsigned int pnf = node2->pn;
          int phase2 = node2->phase;
//        printf("%d %d %d %d\n", m_from_p(pni, ns, wd->n_proton_i, wd->jz_shell), m_from_p(ppi, ns, wd->n_neutron_i, wd->jz_shell), m_from_p(pnf, ns, wd->n_proton_f, wd->jz_shell), m_from_p(ppf, ns, wd->n_neutron_f, wd->jz_shell));
          if (i_op == 0) {
            unsigned int p_hash_i = ppi + wd->n_sds_p_i*(pni % HASH_SIZE);
            unsigned int p_hash_f = ppf + wd->n_sds_p_f*(pnf % HASH_SIZE);
            wh_list *node3 = wd->wh_hash_i[p_hash_i]; // hash corresponds to a_op operators
   
            int index_i = -1;
            int index_f = -1;
   //         int n_hash = 0;
            while (node3 != NULL) {
     //         n_hash++;
              if ((pni == node3->pn) && (ppi == node3->pp)) {
                index_i = node3->index;
                break;
              }
              node3 = node3->next;
            }
  //          printf("n_hash_i: %d\n", n_hash);
            if (index_i < 0) {node2 = node2->next; /*printf("Not found\n"); */continue;}
       //     n_hash = 0;
            node3 = wd->wh_hash_f[p_hash_f];
            while (node3 != NULL) {
         //     n_hash++;
              if ((pnf == node3->pn) && (ppf == node3->pp)) {
                index_f = node3->index;
                break;
              }
              node3 = node3->next;
            }
        //    printf("n_hash_f: %d\n", n_hash);
            if (index_f < 0) {node2 = node2->next; /*printf("Not found\n");*/ continue;}
            total += wd->bc_i[psi_i + wd->n_eig_i*index_i]*wd->bc_f[psi_f + wd->n_eig_f*index_f]*phase1*phase2;
           } else {
            unsigned int p_hash_i = pni + wd->n_sds_p_i*(ppi % HASH_SIZE);
            unsigned int p_hash_f = pnf + wd->n_sds_p_f*(ppf % HASH_SIZE);
            wh_list *node3 = wd->wh_hash_i[p_hash_i]; // hash corresponds to a_op operators
   
            int index_i = -1;
            int index_f = -1;
       //     int n_hash = 0;
            while (node3 != NULL) {
       //       n_hash++;
              if ((pni == node3->pp) && (ppi == node3->pn)) {
                index_i = node3->index;
                break;
              }
              node3 = node3->next;
            }
            if (index_i < 0) {node2 = node2->next; printf("SD not found i\n"); continue;}
      //      printf("N_hash_i: %d\n", n_hash);
      //      n_hash = 0;
            node3 = wd->wh_hash_f[p_hash_f];
            while (node3 != NULL) {
     //         n_hash++;
              if ((pnf == node3->pp) && (ppf == node3->pn)) {
                index_f = node3->index;
                break;
              }
              node3 = node3->next;
            }
     //       printf("N_hash_f: %d\n", n_hash);
            if (index_f < 0) {node2 = node2->next; printf("Not found f: %u\n", p_hash_f); continue;}
            total += wd->bc_i[psi_i + wd->n_eig_i*index_i]*wd->bc_f[psi_f + wd->n_eig_f*index_f]*phase1*phase2;

          }
          node2 = node2->next;
        }
        node1 = node1->next;
      }
    }  
  }
  return total;
}

float trace_a20_nodes(int a, int b, int c, int d, int num_mj, int n_sds_p_int1, int n_sds_n_int1, sd_list** p1_list_i, sd_list** n1_list_i, int* p1_array_f, int* n1_array_f, wfnData* wd, int psi_i, int psi_f) {

  float total = 0.0;
  int ns = wd->n_shells;
  for (int imj = 0; imj < num_mj; imj++) {
    sd_list* node_pi = p1_list_i[imj + num_mj*b]; // Get proton states resulting from p_b |p_i>
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
      sd_list* node_ni = n1_list_i[num_mj - imj + num_mj*d]; // Get neutron states resulting from n_d |n_i>
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
        total += wd->bc_i[psi_i + wd->n_eig_i*index_i]*wd->bc_f[psi_f + wd->n_eig_f*index_f]*phase1*phase2*phase3*phase4;
        node_ni = node_ni->next;
      }
      node_pi = node_pi->next;
    }
  }  
  return total;
}

void build_two_body_jumps_i_and_f(int n_s, int n_p, float mj_min, float mj_max, int num_mj, int n_sds_i, int n_sds_int1, int n_sds_int2, int*a1_array_f, int* a2_array_f, wf_list** a0_list_i, sd_list** a1_list_i, sd_list** a2_list_i, int* jz_shell, int* l_shell) {
/*
  Input(s):
    mj_min_i: 
*/
  for (int j = 1; j <= n_sds_i; j++) {
    int j_min = j_min_from_p(n_s, n_p, j);
    float mj = m_from_p(j, n_s, n_p, jz_shell);
    if ((mj < mj_min) || (mj > mj_max)) {continue;}
    int i_parity = (parity_from_p(j, n_s, n_p, l_shell) + 1)/2;
    int i_mj = mj - mj_min;
    if (a0_list_i[i_parity + 2*i_mj] == NULL) {
      a0_list_i[i_parity + 2*i_mj] = create_wf_node(j, NULL);
    } else {
      wf_append(a0_list_i[i_parity + 2*i_mj], j);
    }
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

void build_two_body_jumps_i(int n_s, int n_p, float mj_min, float mj_max, int num_mj, int n_sds_i, wf_list** a0_list_i, sd_list** a1_list_i, sd_list** a2_list_i, int* jz_shell, int* l_shell) {

  for (int j = 1; j <= n_sds_i; j++) {
    int j_min = j_min_from_p(n_s, n_p, j);
    float mj = m_from_p(j, n_s, n_p, jz_shell);
    if ((mj < mj_min) || (mj > mj_max)) {continue;}
    int i_mj = mj - mj_min;
    int i_parity = (parity_from_p(j, n_s, n_p, l_shell) + 1)/2;
    if (a0_list_i[i_parity + 2*i_mj] == NULL) {
      a0_list_i[i_parity + 2*i_mj] = create_wf_node(j, NULL);
    } else {
      wf_append(a0_list_i[i_parity + 2*i_mj], j);
    }

    for (int b = j_min - 1; b < n_s; b++) {
      int phase1;
      int pn1 = a_op(n_s, n_p, j, b + 1, &phase1, j_min);
      if (pn1 == 0) {continue;}
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
      }
    }
  }

  return;
}

void build_two_body_jumps_f(int n_s, int n_p, float mj_min, float mj_max, int num_mj, int n_sds_f, int n_sds_int1, int n_sds_int2, int* a1_array_f, int* a2_array_f, sd_list** a2_list_f, int* jz_shell, int* l_shell) {
  
  for (int j = 1; j <= n_sds_f; j++) {
    int j_min = j_min_from_p(n_s, n_p, j);
    for (int b = j_min - 1; b < n_s; b++) {
      int phase1;
      int pn1 = a_op(n_s, n_p, j, b + 1, &phase1, j_min);
      if (pn1 == 0) {continue;}
      a1_array_f[(pn1 - 1) + b*n_sds_int1] = j*phase1;
      for (int a = j_min - 1; a < b; a++) {
        int phase2;
        int pn2 = a_op(n_s, n_p - 1, pn1, a + 1, &phase2, j_min);
        if (pn2 == 0) {continue;}
        a2_array_f[(pn2 - 1) + n_sds_int2*(b + a*n_s)] = phase1*phase2*j;
        a2_array_f[(pn2 - 1) + n_sds_int2*(a + b*n_s)] = -phase1*phase2*j;
        float mj = m_from_p(pn2, n_s, n_p - 2, jz_shell);
        if ((mj > mj_max) || (mj < mj_min)) {continue;}
        int i_mj = mj - mj_min;
        int i_parity = (parity_from_p(pn2, n_s, n_p - 2, l_shell) + 1)/2;
        if (a2_list_f[i_parity + 2*(i_mj + num_mj*(b + a*n_s))] == NULL) {
          a2_list_f[i_parity + 2*(i_mj + num_mj*(b + a*n_s))] = create_sd_node(pn2, j, phase1*phase2, NULL);
        } else {
          sd_append(a2_list_f[i_parity + 2*(i_mj + num_mj*(b + a*n_s))], pn2, j, phase1*phase2);
        }
        if (a2_list_f[i_parity + 2*(i_mj + num_mj*(a + b*n_s))] == NULL) {
          a2_list_f[i_parity + 2*(i_mj + num_mj*(a + b*n_s))] = create_sd_node(pn2, j, -phase1*phase2, NULL);
        } else {
          sd_append(a2_list_f[i_parity + 2*(i_mj + num_mj*(a + b*n_s))], pn2, j, -phase1*phase2);
        }
      }
    }
  }

  return;
}
/*
void one_body_density(int j_op, int t_op) {
  // Reads in BIGSTICK basis/wavefunction (.trwfn) files along with
  // orbit definition (.sp) files and constructs the one-body density matrices
  // for each initial and final state eigenfunction
  // Initial and final wave functions must share the same orbit file
  // J_op and T_op are the total spin and isospin of the operator
  wfnData *wd;
  char wfn_file_initial[] = "ne20_basis.trwfn";
  char wfn_file_final[] = "ne20_basis.trwfn";
  char orbit_file[] = "sd.sps";
  wd = read_wfn_data(wfn_file_initial, wfn_file_final, orbit_file);
  FILE* out_file;
  out_file = fopen("ge76_density_1", "w");
  printf("j_op: %d t_op: %d\n", j_op, t_op);

  // Determine the number of intermediate proton SDs
  int* max_state = (int*) malloc(sizeof(int)*(wd->n_proton_i - 1));
  for (int i = 0; i < wd->n_proton_i - 1; i++) {
    max_state[i] = wd->n_shells - (wd->n_proton_i - i - 2);
  }
  int n_sds_p_int = p_step(wd->n_shells, wd->n_proton_i - 1, max_state);
  printf("Num intermediate proton sds: %d\n", n_sds_p_int);
  int n_sds_n_int;
  if (wd->n_neutron_i == 0) {
    n_sds_n_int = 0;
  } else {
    // Determine the number of intermediate neutron SDs
    max_state = realloc(max_state, sizeof(int)*(wd->n_neutron_i - 1));
    for (int i = 0; i < wd->n_neutron_i - 1; i++) {
      max_state[i] = wd->n_shells - (wd->n_neutron_i - i - 2);
    }
    n_sds_n_int = p_step(wd->n_shells, wd->n_neutron_i - 1, max_state);
    free(max_state);
  }
  printf("Num intermediate neutron sds: %d\n", n_sds_n_int);

  sd_list **p_list_i = (sd_list**) malloc(sizeof(sd_list*)*wd->n_shells);
  sd_list **n_list_i = (sd_list**) malloc(sizeof(sd_list*)*wd->n_shells);
  printf("Building initial state proton jumps...\n");
  // For each shell, loop over all proton SDs
  // If c_a | SD > is non-zero, store the resulting SD
  for (int a = 0; a < wd->n_shells; a++) {
    p_list_i[a] = create_sd_node(0,0,1,NULL);
    for (int j = 1; j <= wd->n_sds_p_i; j++) {
      int phase;
      unsigned int pn = a_op(wd->n_shells, wd->n_proton_i, j, a + 1, &phase, 1);
      if (pn == 0) {continue;}
      sd_append(p_list_i[a], j, pn, phase);
    }
  }
  printf("Done.\n");
  printf("Building final state proton jumps...\n");
  int *p_list_f = (int*) calloc(wd->n_shells*n_sds_p_int, sizeof(int));
  for (int a = 0; a < wd->n_shells; a++) {
    for (int j = 1; j <= wd->n_sds_p_f; j++) {
      int phase;
      unsigned int pn = a_op(wd->n_shells, wd->n_proton_f, j, a + 1, &phase, 1);
      if (pn == 0) {continue;}
      p_list_f[(pn - 1) + a*n_sds_p_int] = phase*j;
    }
  }
  printf("Done.\n");

  printf("Building initial state neutron jumps...\n");
  
  for (int a = 0; a < wd->n_shells; a++) {
    n_list_i[a] = create_sd_node(0,0,1,NULL);
    for (int j = 1; j <= wd->n_sds_n_i; j++) {
      int phase;
      int pf = a_op(wd->n_shells, wd->n_neutron_i, j, a + 1, &phase, 1);
      if (pf == 0) {continue;}
      sd_append(n_list_i[a], j, pf, phase);
    }
  }
  printf("Done.\n");
  int *n_list_f = (int*) calloc(wd->n_shells*n_sds_n_int, sizeof(int));
  printf("Building final state neutron jumps...\n");
  
  for (int a = 0; a < wd->n_shells; a++) {
    for (int j = 1; j <= wd->n_sds_n_f; j++) {
      int phase;
      unsigned int pn = a_op(wd->n_shells, wd->n_neutron_f, j, a + 1, &phase, 1);
      if (pn == 0) {continue;}
      n_list_f[(pn - 1) + a*n_sds_n_int] = phase*j;
    }
  }
  printf("Done.\n");
 
  // Loop over initial eigenstates
  for (int psi_i = 1; psi_i < 2; psi_i++) {
    float ji = wd->j_nuc_i[psi_i];
    float ti = wd->t_nuc_i[psi_i];
    float mti = 0.5*(wd->n_proton_i - wd->n_neutron_i);
    // Loop over final eigenstates
    for (int psi_f = 1; psi_f < 2; psi_f++) {
      float cg_j = 0.0;
      float cg_t = 0.0;
      float jf = wd->j_nuc_f[psi_f];
      float tf = wd->t_nuc_f[psi_f];
      float mtf = 0.5*(wd->n_proton_f - wd->n_neutron_f);
      cg_j = clebsch_gordan(j_op, ji, jf, 0, 0, 0);
      if (cg_j == 0.0) {continue;}
      cg_t = clebsch_gordan(t_op, ti, tf, 0, mti, mtf);
      if (cg_t == 0.0) {continue;}
      cg_j *= pow(-1.0, j_op + ji + jf)*sqrt(2*j_op + 1)/sqrt(2*jf + 1);
      cg_t *= pow(-1.0, t_op + ti + tf)*sqrt(2*t_op + 1)/sqrt(2*tf + 1);
      printf("Initial state: # %d J: %g T: %g Final state: # %d J: %g T: %g\n", psi_i + 1, ji, ti, psi_f + 1, jf, tf);
        
      // Loop over final state orbits
      for (int i_orb1 = 0; i_orb1 < wd->n_orbits; i_orb1++) {
        float j1 = wd->j_orb[i_orb1];
        // Loop over initial state orbits
        for (int i_orb2 = 0; i_orb2 < wd->n_orbits; i_orb2++) {
          float j2 = wd->j_orb[i_orb2];
          float total = 0.0;
          if ((j_op > j1 + j2) || (j_op < abs(j1 - j2))) {continue;}
          // Loop over initial state SDs
          for (int b = 0; b < wd->n_shells; b++) {
            if (wd->l_shell[b] != wd->l_orb[i_orb2]) {continue;}
            if (wd->n_shell[b] != wd->n_orb[i_orb2]) {continue;}
            if (wd->j_shell[b]/2.0 != j2) {continue;}
            float mj2 = wd->jz_shell[b]/2.0;
            // Loop over initial state shells
            for (int a = 0; a < wd->n_shells; a++) {
              if (wd->l_shell[a] != wd->l_orb[i_orb1]) {continue;}
              if (wd->n_shell[a] != wd->n_orb[i_orb1]) {continue;}
              if (wd->j_shell[a]/2.0 != j1) {continue;}
              float mj1 = wd->jz_shell[a]/2.0;
              float mt1 = 0.5;
              float mt2 = 0.5;
              float d2 = clebsch_gordan(j1, j2, j_op, mj1, -mj2, 0);
              d2 *= clebsch_gordan(0.5, 0.5, t_op, mt1, -mt2, 0);
              d2 *= pow(-1.0, j2 - mj2 + 0.5 - mt2);
                float d1 = 0.0;
 
              if (d2 != 0.0) {
                sd_list *node = p_list_i[b]->next;
                // Loop over SDs resulting from c_a| SD > 
                while (node != NULL) {
                  // Get the resulting SD
                  int ppn = node->pn;
                  // Check if the resulting SD appears in the 
                  // Corresponding list of SDs resulting from C_b | SD >
                  int ppf = p_list_f[(ppn - 1) + a*n_sds_p_int];
                  if (ppf == 0) {node = node->next; continue;}
                  int ppi = node->pi;
                  int phase1 = node->phase;
                  int phase2 = 1;
                  // Negative SDs indicate negative phase
                  if (ppf < 0) {
                    ppf *= -1;
                    phase2 = -1;
                  }
                  node = node->next;
                  wf_list *node2 = wd->p_hash_i[ppi - 1];
                  wf_list *node3 = wd->p_hash_f[ppf - 1];
                  unsigned int index_i;
                  unsigned int index_f;
                  while (node2 != NULL) {
                    index_i = node2->index;
                    node3 = wd->p_hash_f[ppf - 1];
                    while (node3 != NULL) {
                      index_f = node3->index;
                      if (node2->p == node3->p) {
                        d1 += wd->bc_i[psi_i + wd->n_eig_i*index_i]*wd->bc_f[psi_f + wd->n_eig_f*index_f]*phase1*phase2;
                      }
                      node3 = node3->next;
                    }
                    node2 = node2->next;
                  }
         
                }
                total += d2*d1;
              }
              mt1 = -0.5;
              mt2 = -0.5;
              d2 = clebsch_gordan(j1, j2, j_op, mj1, -mj2, 0);
              d2 *= clebsch_gordan(0.5, 0.5, t_op, mt1, -mt2, 0);
              d2 *= pow(-1.0, j2 - mj2 + 0.5 - mt2);
              if (d2 != 0.0) {

                d1 = 0.0;
                sd_list* node = n_list_i[b]->next;
                while (node != NULL) {
                  int pnn = node->pn;
                  int pnf = n_list_f[(pnn - 1) + a*n_sds_n_int];
                  if (pnf == 0) {node = node->next; continue;}
                  int pni = node->pi;
                  int phase1 = node->phase;
                  int phase2 = 1;
                  if (pnf < 0) {
                    pnf *= -1;
                    phase2 = -1;
                  }
                  node = node->next;
                  wf_list *node2 = wd->n_hash_i[pni - 1];
                  wf_list *node3 = wd->n_hash_f[pnf - 1];
                  unsigned int index_i;
                  unsigned int index_f;
                  while (node2 != NULL) {
                    index_i = node2->index;
                    node3 = wd->n_hash_f[pnf - 1];
                    while (node3 != NULL) {
                      index_f = node3->index;
                      if (node2->p == node3->p) {
                        d1 += wd->bc_i[psi_i + wd->n_eig_i*index_i]*wd->bc_f[psi_f + wd->n_eig_f*index_f]*phase1*phase2;
                      }
                      node3 = node3->next;
                    }
                    node2 = node2->next;
                  }
         
                }
                total += d2*d1;
              }
              mt1 = 0.5;
              mt2 = -0.5;
              d2 = clebsch_gordan(j1, j2, j_op, mj1, -mj2, 0);
              d2 *= clebsch_gordan(0.5, 0.5, t_op, mt1, -mt2, 0);
              d2 *= pow(-1.0, j2 - mj2 + 0.5 - mt2);
              if (d2 != 0.0) {

                d1 = 0.0;
                sd_list* node = p_list_i[b]->next;
                while (node != NULL) {
                  int ppf = node->pn;
                  int ppi = node->pi; 
                  int phase1 = node->phase;
                  // Get list of n_f associated to p_f
                  wf_list *node2 = wd->p_hash_f[ppf - 1];
                  unsigned int index_i;
                  unsigned int index_f;
                  // Loop over n_f 
                  while (node2 != NULL) {
                    index_i = node2->index;
                    unsigned int pnf = node2->p;
                    wf_list *node3 = wd->p_hash_i[ppi - 1];
                    while (node3 != NULL) {
                      unsigned int pni = node3->p;
                      int phase2 = 1;
                      if (n_list_f[(pnf - 1) + a*n_sds_n_int] == pni) {
                        index_f = node3->index;
                        d1 += wd->bc_i[psi_i + wd->n_eig_i*index_i]*wd->bc_f[psi_f + wd->n_eig_f*index_f]*phase1*phase2;
                      } else if (n_list_f[(pnf - 1) + a*n_sds_n_int] == -pni) {
                        phase2 = -1;
                        index_f = node3->index;
                        d1 += wd->bc_i[psi_i + wd->n_eig_i*index_i]*wd->bc_f[psi_f + wd->n_eig_f*index_f]*phase1*phase2;
                      }

                      node3 = node3->next;
                    }
                    node2 = node2->next;
                  }
                }  
              }
              total += d1*d2;
            }
          }
          total /= cg_j*cg_t;
          if (total != 0.0) {
            printf("%d %d %g %d %d %g %g\n", wd->n_orb[i_orb1], wd->l_orb[i_orb1], wd->j_orb[i_orb1], wd->n_orb[i_orb2], wd->l_orb[i_orb2], wd->j_orb[i_orb2], total);
          //  printf("%d %d %g\n", i_orb1 + 1, i_orb2 + 1, total);
          }
        }
      }
    }
  }          
  return;
  fclose(out_file);
}
*/

