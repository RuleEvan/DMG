#include "file_io.h"

int VERBOSE = I_VERBOSE;

speedParams* read_parameter_file(char *param_file) {
  FILE *in_file;
  speedParams *sp = malloc(sizeof(*sp));
  in_file = fopen(param_file, "r");
  sp->initial_file_base = malloc(sizeof(char)*100);
  fscanf(in_file, "%s\n", sp->initial_file_base);
  sp->final_file_base = malloc(sizeof(char)*100);
  fscanf(in_file, "%s\n", sp->final_file_base);
  sp->out_file_base = malloc(sizeof(char)*100);
  fscanf(in_file, "%s\n", sp->out_file_base);
  fscanf(in_file, "%d\n", &sp->n_body);
  fscanf(in_file, "%d,%d\n", &sp->j_op, &sp->t_op);
  if ((sp->n_body != 1) && (sp->n_body != 2)) {printf("Incorrect request for %d body operator, please type only 1 or 2 here.\n", sp->n_body); exit(0);}
  int eig_i, eig_f;
  sp->n_trans = 0;
  sp->transition_list = NULL;
  while (fscanf(in_file, "%d,%d\n", &eig_i, &eig_f) == 2) {
    (sp->n_trans)++;
    if (sp->transition_list == NULL) {
      sp->transition_list = create_eigen_node(eig_i, eig_f, NULL);
    } else {
      eigen_append(sp->transition_list, eig_i, eig_f);
    }
  }
  if (VERBOSE) {printf("Read in %d transitions\n", sp->n_trans);}
  return sp;
}

wfnData* read_binary_wfn_data(char *wfn_file_initial, char *wfn_file_final, char *basis_file_initial, char *basis_file_final) {
  wfnData *wd = malloc(sizeof(*wd));
  FILE *in_file;
  // Read in initial wavefunction data
  if (VERBOSE) {printf("\n");
  printf("Opening initial state wfn file\n");}
  in_file = fopen(wfn_file_initial, "rb");
  int junk;
  fread(&junk, sizeof(int), 1, in_file);
  fread(&junk, sizeof(int), 1, in_file);
  int vec_offset;
  char buffer[100];
  fread(&vec_offset, sizeof(int), 1, in_file);
  fread(&junk, sizeof(int), 1, in_file);
  fread(&junk, sizeof(int), 1, in_file);
  fread(&junk, sizeof(int), 1, in_file);
  fread(&junk, sizeof(int), 1, in_file);
  fread(&junk, sizeof(int), 1, in_file);
  fread(&wd->n_proton_i, sizeof(int), 1, in_file);
  fread(&wd->n_neutron_i, sizeof(int), 1, in_file);
  wd->n_data = wd->n_proton_i + wd->n_neutron_i;
  if (VERBOSE) {printf("Initial state contains %d protons and %d neutrons\n", wd->n_proton_i, wd->n_neutron_i);}

  fread(&junk, sizeof(int), 1, in_file);
//  fread(&junk, sizeof(int), 1, in_file);
  int n_proton_orbs, n_neutron_orbs;
  fread(&n_proton_orbs, sizeof(int), 1, in_file);
  fread(&n_neutron_orbs, sizeof(int), 1, in_file);
  wd->n_orbits = n_proton_orbs;
  if (VERBOSE) {printf("Model space contains %d orbitals\n", n_neutron_orbs);}
  wd->n_orb = (int*) malloc(sizeof(int)*wd->n_orbits);
  wd->l_orb = (int*) malloc(sizeof(int)*wd->n_orbits);
  wd->w_orb = (int*) malloc(sizeof(int)*wd->n_orbits);
  wd->j_orb = (float*) malloc(sizeof(float)*wd->n_orbits);
  for (int i = 0; i < n_proton_orbs; i++) {
    int j_orb, pi_orb, w_orb;
    fread(&wd->n_orb[i], sizeof(int), 1, in_file);
    fread(&j_orb, sizeof(int), 1, in_file);
    wd->j_orb[i] = j_orb/2.0;
    fread(&wd->l_orb[i], sizeof(int), 1, in_file);
    fread(&pi_orb, sizeof(int), 1, in_file);
    fread(&wd->w_orb[i], sizeof(int), 1, in_file);
  }
  for (int i = 0; i < n_neutron_orbs; i++) {
    int n_orb, j_orb, l_orb, pi_orb, w_orb;
    fread(&n_orb, sizeof(int), 1, in_file);
    fread(&j_orb, sizeof(int), 1, in_file);
    fread(&l_orb, sizeof(int), 1, in_file);
    fread(&pi_orb, sizeof(int), 1, in_file);
    fread(&w_orb, sizeof(int), 1, in_file);
  }
  int n_sps_p, n_sps_n;
  fread(&n_sps_p, sizeof(int), 1, in_file);
  fread(&n_sps_n, sizeof(int), 1, in_file);
  wd->n_shells = n_sps_p;
  if (VERBOSE) {printf("Number of single particle states: %d\n", wd->n_shells);}
  wd->n_shell = (int*) malloc(sizeof(int)*wd->n_shells);
  wd->j_shell = (int*) malloc(sizeof(int)*wd->n_shells);
  wd->l_shell = (int*) malloc(sizeof(int)*wd->n_shells);
  wd->jz_shell = (int*) malloc(sizeof(int)*wd->n_shells);
  wd->w_shell = (int*) malloc(sizeof(int)*wd->n_shells);
  wd->n_sds_p_i = get_num_sds(wd->n_shells, wd->n_proton_i);
  wd->n_sds_n_i = get_num_sds(wd->n_shells, wd->n_neutron_i);
  if (VERBOSE) {printf("Initial proton SDs: %d Initial neutron SDs: %d\n", wd->n_sds_p_i, wd->n_sds_n_i);}
  for (int i = 0; i < n_sps_p; i++) {
    int n_shell, l_shell, j_shell, jz_shell, w_shell, pi_shell, id_shell, gid_shell;
    fread(&n_shell, sizeof(int), 1, in_file);
    fread(&j_shell, sizeof(int), 1, in_file);
    fread(&jz_shell, sizeof(int), 1, in_file);
    fread(&l_shell, sizeof(int), 1, in_file);
    fread(&w_shell, sizeof(int), 1, in_file);
    fread(&pi_shell, sizeof(int), 1, in_file);
    fread(&id_shell, sizeof(int), 1, in_file);
    fread(&gid_shell, sizeof(int), 1, in_file);
  }
  for (int i = 0; i < n_sps_n; i++) {
    int n_shell, j_shell, jz_shell, l_shell, w_shell, pi_shell, id_shell, gid_shell;

    fread(&n_shell, sizeof(int), 1, in_file);
    fread(&j_shell, sizeof(int), 1, in_file);
    fread(&jz_shell, sizeof(int), 1, in_file);
    fread(&l_shell, sizeof(int), 1, in_file);
    fread(&w_shell, sizeof(int), 1, in_file);
    fread(&pi_shell, sizeof(int), 1, in_file);
    fread(&id_shell, sizeof(int), 1, in_file);
    fread(&gid_shell, sizeof(int), 1, in_file);
  }


  fread(&wd->jz_i, sizeof(int), 1, in_file);
  fread(&wd->parity_i, sizeof(int), 1, in_file);
  fread(&wd->w_max_i, sizeof(int), 1, in_file);
  fread(&junk, sizeof(int), 1, in_file);
  fread(&junk, sizeof(int), 1, in_file);
  fread(&junk, sizeof(int), 1, in_file);
  fread(&wd->n_states_i, sizeof(long long int), 1, in_file);
  if (VERBOSE) {printf("The model space has %d shells for a total basis size of %d\n", wd->n_shells, wd->n_states_i);}
  fread(&wd->n_eig_i, sizeof(int), 1, in_file);
  if (VERBOSE) {printf("Initial state file contains %d eigenstates\n", wd->n_eig_i);}
  if (wd->parity_i == '+') {if (VERBOSE) {printf("The initial state basis is positive parity\n");}} 
  else if (wd->parity_i == '-') {if (VERBOSE) {printf("The initial state basis is negative parity\n");}} 
  else if (wd->parity_i == '0') {if (VERBOSE) {printf("The initial state basis contains both positive and negative parity states\n");}} 
  else {printf("Error in initial state parity: %c\n", wd->parity_i); exit(0);}
  if (wd->w_max_i > 0) {
    if (VERBOSE) {
      printf("The max initial state excitation is w_max = %d\n", wd->w_max_i);
    }
  } else if (wd->w_max_i < 0) {
    printf("Error: Max excitation w_max < 0. Please define a truncation scheme where w_max is positive\n"); 
    exit(0);
  }

  wd->e_nuc_i =  (float*) malloc(sizeof(float)*wd->n_eig_i);
  wd->j_nuc_i =  (float*) malloc(sizeof(float)*wd->n_eig_i);
  wd->t_nuc_i =  (float*) malloc(sizeof(float)*wd->n_eig_i);
  unsigned int wave_size = wd->n_states_i*wd->n_eig_i;
  wd->bc_i = (float*) calloc(wave_size, sizeof(float));
  
  if (VERBOSE) {printf("Reading in initial state wavefunction coefficients\n");}
  for (unsigned int i = 0; i < wd->n_eig_i; i++) {
    fread(&junk, sizeof(int), 1, in_file);
    int vec_index;
    double total = 0.0;
    float j_nuc, t_nuc;
    fread(&vec_index, sizeof(int), 1, in_file);
    fread(&wd->e_nuc_i[i], sizeof(float), 1, in_file);
    fread(&j_nuc, sizeof(float), 1, in_file);
    fread(&t_nuc, sizeof(float), 1, in_file);
    t_nuc = -0.5 + 0.5*sqrt(1 + 4.0*t_nuc);
    int ij_nuc = round(2*j_nuc);
    int it_nuc = round(2*t_nuc);
    wd->j_nuc_i[i] = ij_nuc/2.0;
    wd->t_nuc_i[i] = it_nuc/2.0;
    for (unsigned int j = 0; j < wd->n_states_i; j++) {
      unsigned int index = i + wd->n_eig_i*j;
      //fread(&wd->bc_i[i + wd->n_eig_i*j], sizeof(float), 1, in_file);
     // total += pow(wd->bc_i[i + wd->n_eig_i*j], 2);
      fread(&wd->bc_i[index], sizeof(float), 1, in_file);
      total += pow(wd->bc_i[index], 2);

    }
    if (fabs(total -1.0) > pow(10, -6)) {printf("State %d not normalized: norm = %g\n", i, total); exit(0);}
  }
  fclose(in_file);
  if (VERBOSE) {printf("Done.\n");
  printf("Reading in initial state basis\n");}
  in_file = fopen(basis_file_initial, "rb");

  fread(&junk, sizeof(int), 1, in_file);
  fread(&junk, sizeof(int), 1, in_file);
  fread(&junk, sizeof(int), 1, in_file);
  fread(&vec_offset, sizeof(int), 1, in_file);
  fseek(in_file, vec_offset + 16, SEEK_SET);
  for (int i = 0; i < wd->n_shells; i++) {
    fread(&junk, sizeof(int), 1, in_file);
    fread(&wd->n_shell[i], sizeof(int), 1, in_file);
    fread(&wd->l_shell[i], sizeof(int), 1, in_file);
    fread(&wd->j_shell[i], sizeof(int), 1, in_file);
    fread(&wd->jz_shell[i], sizeof(int), 1, in_file);
    fread(&wd->w_shell[i], sizeof(int), 1, in_file);
  }
  for (int i = 0; i < wd->n_shells; i++) {
    int n_shell, j_shell, l_shell, jz_shell, w_shell;
    fread(&junk, sizeof(int), 1, in_file);
    fread(&n_shell, sizeof(int), 1, in_file);
    fread(&l_shell, sizeof(int), 1, in_file);
    fread(&j_shell, sizeof(int), 1, in_file);
    fread(&jz_shell, sizeof(int), 1, in_file);
    fread(&w_shell, sizeof(int), 1, in_file);
  }
  
  wd->wh_hash_i = (wh_list**) calloc(wd->n_sds_p_i*HASH_SIZE, sizeof(wh_list*));
  unsigned int pp0 = n_choose_k(wd->n_shells, wd->n_proton_i);
  unsigned int pn0 = n_choose_k(wd->n_shells, wd->n_neutron_i);
  
  for (int i = 0; i < wd->n_states_i; i++) {
    int in = 0;
    int ip = 0;
    unsigned int pp = pp0;
    unsigned int pn = pn0;
    for (int j = 0; j < wd->n_data; j++) {
      int i_state;
      fread(&i_state, sizeof(int), 1, in_file);
      if (i_state <= wd->n_shells) {
        pp -= n_choose_k(wd->n_shells - i_state, wd->n_proton_i - ip);
        ip++;
      } else {
        pn -= n_choose_k(2*(wd->n_shells) - i_state, wd->n_neutron_i - in);
        in++;
      }
    }
    //printf("%u, %u\n", pp, pn);
    unsigned int p_hash = pp + wd->n_sds_p_i*(pn % HASH_SIZE);
    
    if (wd->wh_hash_i[p_hash] == NULL) {
      wd->wh_hash_i[p_hash] = create_wh_node(pp, pn, i, NULL);
    } else {
      wh_append(wd->wh_hash_i[p_hash], pp, pn, i);
    }
  }
  if (VERBOSE) {printf("Done.\n");}
  fclose(in_file);
  if (VERBOSE) {printf("\n");}
  if (strcmp(wfn_file_initial, wfn_file_final) == 0) {
    if (VERBOSE) {printf("Initial and final bases are identical.\n");}
    wd->same_basis = 1;
    wd->n_proton_f = wd->n_proton_i;
    wd->n_neutron_f = wd->n_neutron_i;
    wd->n_states_f = wd->n_states_i;
    wd->n_sds_p_f = wd->n_sds_p_i;
    wd->n_sds_n_f = wd->n_sds_n_i;
    wd->jz_f = wd->jz_i;
    wd->n_eig_f = wd->n_eig_i;    
    wd->e_nuc_f =  (float*) malloc(sizeof(float)*wd->n_eig_f);
    wd->j_nuc_f =  (float*) malloc(sizeof(float)*wd->n_eig_f);
    wd->t_nuc_f =  (float*) malloc(sizeof(float)*wd->n_eig_f);
    wd->e_nuc_f = wd->e_nuc_i;
    wd->j_nuc_f = wd->j_nuc_i;
    wd->t_nuc_f = wd->t_nuc_i;
    wd->bc_f = wd->bc_i;
    wd->wh_hash_f = wd->wh_hash_i;

  } else {
    if (VERBOSE) {printf("Initial and final bases are different.\n");}
    wd->same_basis = 0;
    in_file = fopen(wfn_file_final, "rb");
    fread(&junk, sizeof(int), 1, in_file);
    fread(&junk, sizeof(int), 1, in_file);
    fread(&vec_offset, sizeof(int), 1, in_file);
    fread(&junk, sizeof(int), 1, in_file);
    fread(&junk, sizeof(int), 1, in_file);
    fread(&junk, sizeof(int), 1, in_file);
    fread(&junk, sizeof(int), 1, in_file);
    fread(&junk, sizeof(int), 1, in_file);
    fread(&wd->n_proton_f, sizeof(int), 1, in_file);
    fread(&wd->n_neutron_f, sizeof(int), 1, in_file);
    if (VERBOSE) {printf("Final state contains %d protons and %d neutrons\n", wd->n_proton_f, wd->n_neutron_f);}
    wd->n_sds_p_f = get_num_sds(wd->n_shells, wd->n_proton_f);
    wd->n_sds_n_f = get_num_sds(wd->n_shells, wd->n_neutron_f);
    if (VERBOSE) {printf("Final proton SDs: %d Final neutron SDs: %d\n", wd->n_sds_p_f, wd->n_sds_n_f);}
    fread(&junk, sizeof(int), 1, in_file);
    fread(&n_proton_orbs, sizeof(int), 1, in_file);
    fread(&n_neutron_orbs, sizeof(int), 1, in_file);
    for (int i = 0; i < n_proton_orbs; i++) {
      int j_orb, pi_orb, w_orb;
      fread(&wd->n_orb[i], sizeof(int), 1, in_file);
      fread(&j_orb, sizeof(int), 1, in_file);
      wd->j_orb[i] = j_orb/2.0;
      fread(&wd->l_orb[i], sizeof(int), 1, in_file);
      fread(&pi_orb, sizeof(int), 1, in_file);
      fread(&w_orb, sizeof(int), 1, in_file);
    }
    for (int i = 0; i < n_neutron_orbs; i++) {
      int n_orb, j_orb, l_orb, pi_orb, w_orb;
      fread(&n_orb, sizeof(int), 1, in_file);
      fread(&j_orb, sizeof(int), 1, in_file);
      fread(&l_orb, sizeof(int), 1, in_file);
      fread(&pi_orb, sizeof(int), 1, in_file);
      fread(&w_orb, sizeof(int), 1, in_file);
    }
    fread(&n_sps_p, sizeof(int), 1, in_file);
    fread(&n_sps_n, sizeof(int), 1, in_file);
    for (int i = 0; i < n_sps_p; i++) {
    int n_shell, j_shell, jz_shell, l_shell, w_shell, pi_shell, id_shell, gid_shell;

      fread(&n_shell, sizeof(int), 1, in_file);
      fread(&j_shell, sizeof(int), 1, in_file);
      fread(&jz_shell, sizeof(int), 1, in_file);
      fread(&l_shell, sizeof(int), 1, in_file);
      fread(&w_shell, sizeof(int), 1, in_file);
      fread(&pi_shell, sizeof(int), 1, in_file);
      fread(&id_shell, sizeof(int), 1, in_file);
      fread(&gid_shell, sizeof(int), 1, in_file);
    }
    for (int i = 0; i < n_sps_n; i++) {
      int n_shell, j_shell, jz_shell, l_shell, w_shell, pi_shell, id_shell, gid_shell;

      fread(&n_shell, sizeof(int), 1, in_file);
      fread(&j_shell, sizeof(int), 1, in_file);
      fread(&jz_shell, sizeof(int), 1, in_file);
      fread(&l_shell, sizeof(int), 1, in_file);
      fread(&w_shell, sizeof(int), 1, in_file);
      fread(&pi_shell, sizeof(int), 1, in_file);
      fread(&id_shell, sizeof(int), 1, in_file);
      fread(&gid_shell, sizeof(int), 1, in_file);
    }

    fread(&wd->jz_f, sizeof(int), 1, in_file);
    fread(&wd->parity_f, sizeof(int), 1, in_file);
    fread(&wd->w_max_f, sizeof(int), 1, in_file);
    fread(&junk, sizeof(int), 1, in_file);
    fread(&junk, sizeof(int), 1, in_file);
    fread(&junk, sizeof(int), 1, in_file);
    fread(&wd->n_states_f, sizeof(long long int), 1, in_file);
    if (VERBOSE) {printf("The model space has %d shells for a total basis size of %d\n", wd->n_shells, wd->n_states_f);}
    fread(&wd->n_eig_f, sizeof(int), 1, in_file);
    if (VERBOSE) {printf("Final state file contains %d eigenstates\n", wd->n_eig_f);}
 
    if (wd->parity_f == '+') {if (VERBOSE) {printf("The final state basis is positive parity\n");}} 
    else if (wd->parity_f == '-') {if (VERBOSE) {printf("The final state basis is negative parity\n");}} 
    else if (wd->parity_f == '0') {if (VERBOSE) {printf("The final state basis contains both positive and negative parity states\n");}} 
    else {printf("Error in final state parity: %c\n", wd->parity_f); exit(0);}
    if (wd->w_max_f > 0) {
      if (VERBOSE) {printf("The max final state excitation is w_max = %d\n", wd->w_max_f);}
    } else if (wd->w_max_f < 0) {
      printf("Error: Max excitation w_max < 0. Please define a truncation scheme where w_max > 0.\n");
      exit(0);
    }

    wd->e_nuc_f = (float*) malloc(sizeof(float)*wd->n_eig_f);
    wd->j_nuc_f = (float*) malloc(sizeof(float)*wd->n_eig_f);
    wd->t_nuc_f = (float*) malloc(sizeof(float)*wd->n_eig_f);
    wave_size = wd->n_states_f*wd->n_eig_f;
    wd->bc_f = (float*) calloc(wave_size, sizeof(float));

    //wd->bc_f = malloc(sizeof(float)*wd->n_states_f*wd->n_eig_f);

    if (VERBOSE) {printf("Reading in final state wavefunction coefficients\n");}
    for (unsigned int i = 0; i < wd->n_eig_f; i++) {
      fread(&junk, sizeof(int), 1, in_file);
      int vec_index;
      float j_nuc, t_nuc;
      double total = 0.0;
      fread(&vec_index, sizeof(int), 1, in_file);
      fread(&wd->e_nuc_f[i], sizeof(float), 1, in_file);
      fread(&j_nuc, sizeof(float), 1, in_file);
      fread(&t_nuc, sizeof(float), 1, in_file);
      t_nuc = -0.5 + 0.5*sqrt(1 + 4.0*t_nuc);
      int ij_nuc = round(2*j_nuc);
      int it_nuc = round(2*t_nuc);
      wd->j_nuc_f[i] = ij_nuc/2.0;
      wd->t_nuc_f[i] = it_nuc/2.0;
      for (unsigned int j = 0; j < wd->n_states_f; j++) {
        fread(&wd->bc_f[i + wd->n_eig_f*j], sizeof(float), 1, in_file);
        total += pow(wd->bc_f[i + wd->n_eig_f*j], 2.0);
      }
      if (fabs(total - 1.0) > pow(10, -6)) {printf("State %d is not normalized: norm = %g\n", i, total); exit(0);}
    }
    fclose(in_file);
    if (VERBOSE) {printf("Done.\n");
    printf("Reading in final state basis\n");}
    in_file = fopen(basis_file_final, "rb");
    fread(&junk, sizeof(int), 1, in_file);
    fread(&junk, sizeof(int), 1, in_file);
    fread(&junk, sizeof(int), 1, in_file);
    fread(&vec_offset, sizeof(int), 1, in_file);
    fseek(in_file, vec_offset + 16, SEEK_SET);
    for (int i = 0; i < 2*wd->n_shells; i++) {
      int n_shell, j_shell, l_shell, jz_shell, w_shell;
      fread(&junk, sizeof(int), 1, in_file);
      fread(&n_shell, sizeof(int), 1, in_file);
      fread(&l_shell, sizeof(int), 1, in_file);
      fread(&j_shell, sizeof(int), 1, in_file);
      fread(&jz_shell, sizeof(int), 1, in_file);
      fread(&w_shell, sizeof(int), 1, in_file);
      //printf("%d, %d, %d, %d, %d\n", n_shell, l_shell, j_shell, jz_shell, w_shell);
    }
    wd->wh_hash_f = (wh_list**) calloc(wd->n_sds_p_f*HASH_SIZE, sizeof(wh_list*));
    int pp0 = n_choose_k(wd->n_shells, wd->n_proton_f);
    int pn0 = n_choose_k(wd->n_shells, wd->n_neutron_f);
    for (int i = 0; i < wd->n_states_f; i++) {
      int in = 0;
      int ip = 0;
      unsigned int pp = pp0;
      unsigned int pn = pn0;
      for (int j = 0; j < wd->n_data; j++) {
        int i_state;
        fread(&i_state, sizeof(int), 1, in_file);
        if (i_state <= wd->n_shells) {
          pp -= n_choose_k(wd->n_shells - i_state, wd->n_proton_f - ip);
          ip++;
        } else {
          pn -= n_choose_k(2*(wd->n_shells) - i_state, wd->n_neutron_f - in);
          in++;
        }
      }
      unsigned int p_hash = pp + wd->n_sds_p_f*(pn % HASH_SIZE);
      if (wd->wh_hash_f[p_hash] == NULL) {
        wd->wh_hash_f[p_hash] = create_wh_node(pp, pn, i, NULL);
      } else {
        wh_append(wd->wh_hash_f[p_hash], pp, pn, i);
      }
    }
    if (VERBOSE) {printf("Done.\n");}
    fclose(in_file);

  }

  return wd;
}

sd_list* create_sd_node(unsigned int pi, unsigned int pn, int phase, sd_list* next) {
  sd_list* new_node = (sd_list*)malloc(sizeof(sd_list));
  if (new_node == NULL) {
    printf("Error creating node\n");
    exit(0);
  }
  new_node->pi = pi;
  new_node->pn = pn;
  new_node->phase = phase;
  new_node->next = next;

  return new_node;
}


sd_list* sd_append(sd_list* head, unsigned int pi, unsigned int pn, int phase) {
  if (head->next == NULL) {
    sd_list* new_node = create_sd_node(pi, pn, phase, NULL);
    head->next = new_node;
  } else {
    sd_list* next = head->next;
    sd_list* new_node = create_sd_node(pi, pn, phase, next);
    head->next = new_node;
  }

  
  return head;
}

sde_list* create_sde_node(unsigned int pi, unsigned int pn, int phase, int n_quanta, sde_list* next) {
  sde_list* new_node = (sde_list*)malloc(sizeof(sde_list));
  if (new_node == NULL) {
    printf("Error creating node\n");
    exit(0);
  }
  new_node->pi = pi;
  new_node->pn = pn;
  new_node->phase = phase;
  new_node->n_quanta = n_quanta;
  new_node->next = next;

  return new_node;
}

sde_list* sde_append(sde_list* head, unsigned int pi, unsigned int pn, int phase, int n_quanta) {
  if (head->next == NULL) {
    sde_list* new_node = create_sde_node(pi, pn, phase, n_quanta, NULL);
    head->next = new_node;
  } else {
    sde_list* next = head->next;
    sde_list* new_node = create_sde_node(pi, pn, phase, n_quanta, next);
    head->next = new_node;
  }

  
  return head;
}


wf_list* create_wf_node(unsigned int p, wf_list* next) {
  wf_list* new_node = (wf_list*)malloc(sizeof(wf_list));
  if (new_node == NULL) {
    printf("Error creating node\n");
    exit(0);
  }
  new_node->p = p;
  new_node->next = next;
  
  return new_node;
}

wf_list* wf_append(wf_list* head, unsigned int p) {
  if (head->next == NULL) {
    wf_list* new_node = create_wf_node(p, NULL);
    head->next = new_node;
  } else {
    wf_list* next = head->next;
    wf_list* new_node = create_wf_node(p, next);
    head->next = new_node;
  }

  return head;
}

wfe_list* create_wfe_node(unsigned int p, int n_quanta, wfe_list* next) {
  wfe_list* new_node = (wfe_list*)malloc(sizeof(wfe_list));
  if (new_node == NULL) {
    printf("Error creating node\n");
    exit(0);
  }
  new_node->p = p;
  new_node->n_quanta = n_quanta;
  new_node->next = next;
  
  return new_node;
}

wfe_list* wfe_append(wfe_list* head, unsigned int p, int n_quanta) {
  if (head->next == NULL) {
    wfe_list* new_node = create_wfe_node(p, n_quanta, NULL);
    head->next = new_node;
  } else {
    wfe_list* next = head->next;
    wfe_list* new_node = create_wfe_node(p, n_quanta, next);
    head->next = new_node;
  }

  return head;
}


wh_list* create_wh_node(unsigned int pp, unsigned int pn, unsigned int index, wh_list* next) {
  wh_list* new_node = (wh_list*)malloc(sizeof(wh_list));
  if (new_node == NULL) {
    printf("Error creating node\n");
    exit(0);
  }
  new_node->pp = pp;
  new_node->pn = pn;
  new_node->index = index;
  new_node->next = next;
  
  return new_node;
}

wh_list* wh_append(wh_list* head, unsigned int pp, unsigned int pn, unsigned int index) {
  if (head->next == NULL) {
    wh_list* new_node = create_wh_node(pp, pn, index, NULL);
    head->next = new_node;
  } else {
    wh_list* next = head->next;
    wh_list* new_node = create_wh_node(pp, pn, index, next);
    head->next = new_node;
  }

  return head;
}

eigen_list* create_eigen_node(int eig_i, int eig_f, eigen_list* next) {
  eigen_list* new_node = (eigen_list*)malloc(sizeof(eigen_list));
  if (new_node == NULL) {
    printf("Error creating node\n");
    exit(0);
  }
  new_node->eig_i = eig_i;
  new_node->eig_f = eig_f;
  new_node->next = next;
  
  return new_node;
}

eigen_list* eigen_append(eigen_list* head, int eig_i, int eig_f) {
  if (head->next == NULL) {
    eigen_list* new_node = create_eigen_node(eig_i, eig_f, NULL);
    head->next = new_node;
  } else {
    head = eigen_append(head->next, eig_i, eig_f);
  }

  return head;
}

wfnData* read_binary_wfn_strength(char *wfn_file_initial, char *wfn_file_final, char *basis_file_initial, char *basis_file_final) {
  wfnData *wd = malloc(sizeof(*wd));
  FILE *in_file;
  // Read in initial wavefunction data
  if (VERBOSE) {printf("\n");
  printf("Opening initial state wfn file\n");}
  in_file = fopen(wfn_file_initial, "rb");
  int junk;
  fread(&junk, sizeof(int), 1, in_file);
  fread(&junk, sizeof(int), 1, in_file);
  int vec_offset;
  char buffer[100];
  fread(&vec_offset, sizeof(int), 1, in_file);
  fread(&junk, sizeof(int), 1, in_file);
  fread(&junk, sizeof(int), 1, in_file);
  fread(&junk, sizeof(int), 1, in_file);
  fread(&junk, sizeof(int), 1, in_file);
  fread(&junk, sizeof(int), 1, in_file);
  fread(&wd->n_proton_i, sizeof(int), 1, in_file);
  fread(&wd->n_neutron_i, sizeof(int), 1, in_file);
  wd->n_data = wd->n_proton_i + wd->n_neutron_i;
  if (VERBOSE) {printf("Initial state contains %d protons and %d neutrons\n", wd->n_proton_i, wd->n_neutron_i);}

  fread(&junk, sizeof(int), 1, in_file);
//  fread(&junk, sizeof(int), 1, in_file);
  int n_proton_orbs, n_neutron_orbs;
  fread(&n_proton_orbs, sizeof(int), 1, in_file);
  fread(&n_neutron_orbs, sizeof(int), 1, in_file);
  wd->n_orbits = n_proton_orbs;
  if (VERBOSE) {printf("Model space contains %d orbitals\n", n_neutron_orbs);}
  wd->n_orb = (int*) malloc(sizeof(int)*wd->n_orbits);
  wd->l_orb = (int*) malloc(sizeof(int)*wd->n_orbits);
  wd->w_orb = (int*) malloc(sizeof(int)*wd->n_orbits);
  wd->j_orb = (float*) malloc(sizeof(float)*wd->n_orbits);
  for (int i = 0; i < n_proton_orbs; i++) {
    int j_orb, pi_orb, w_orb;
    fread(&wd->n_orb[i], sizeof(int), 1, in_file);
    fread(&j_orb, sizeof(int), 1, in_file);
    wd->j_orb[i] = j_orb/2.0;
    fread(&wd->l_orb[i], sizeof(int), 1, in_file);
    fread(&pi_orb, sizeof(int), 1, in_file);
    fread(&wd->w_orb[i], sizeof(int), 1, in_file);
  }
  for (int i = 0; i < n_neutron_orbs; i++) {
    int n_orb, j_orb, l_orb, pi_orb, w_orb;
    fread(&n_orb, sizeof(int), 1, in_file);
    fread(&j_orb, sizeof(int), 1, in_file);
    fread(&l_orb, sizeof(int), 1, in_file);
    fread(&pi_orb, sizeof(int), 1, in_file);
    fread(&w_orb, sizeof(int), 1, in_file);
  }
  int n_sps_p, n_sps_n;
  fread(&n_sps_p, sizeof(int), 1, in_file);
  fread(&n_sps_n, sizeof(int), 1, in_file);
  wd->n_shells = n_sps_p;
  if (VERBOSE) {printf("Number of single particle states: %d\n", wd->n_shells);}
  wd->n_shell = (int*) malloc(sizeof(int)*wd->n_shells);
  wd->j_shell = (int*) malloc(sizeof(int)*wd->n_shells);
  wd->l_shell = (int*) malloc(sizeof(int)*wd->n_shells);
  wd->jz_shell = (int*) malloc(sizeof(int)*wd->n_shells);
  wd->w_shell = (int*) malloc(sizeof(int)*wd->n_shells);
  wd->n_sds_p_i = get_num_sds(wd->n_shells, wd->n_proton_i);
  wd->n_sds_n_i = get_num_sds(wd->n_shells, wd->n_neutron_i);
  if (VERBOSE) {printf("Initial proton SDs: %d Initial neutron SDs: %d\n", wd->n_sds_p_i, wd->n_sds_n_i);}
  for (int i = 0; i < n_sps_p; i++) {
    int n_shell, l_shell, j_shell, jz_shell, w_shell, pi_shell, id_shell, gid_shell;
    fread(&n_shell, sizeof(int), 1, in_file);
    fread(&j_shell, sizeof(int), 1, in_file);
    fread(&jz_shell, sizeof(int), 1, in_file);
    fread(&l_shell, sizeof(int), 1, in_file);
    fread(&w_shell, sizeof(int), 1, in_file);
    fread(&pi_shell, sizeof(int), 1, in_file);
    fread(&id_shell, sizeof(int), 1, in_file);
    fread(&gid_shell, sizeof(int), 1, in_file);
  }
  for (int i = 0; i < n_sps_n; i++) {
    int n_shell, j_shell, jz_shell, l_shell, w_shell, pi_shell, id_shell, gid_shell;

    fread(&n_shell, sizeof(int), 1, in_file);
    fread(&j_shell, sizeof(int), 1, in_file);
    fread(&jz_shell, sizeof(int), 1, in_file);
    fread(&l_shell, sizeof(int), 1, in_file);
    fread(&w_shell, sizeof(int), 1, in_file);
    fread(&pi_shell, sizeof(int), 1, in_file);
    fread(&id_shell, sizeof(int), 1, in_file);
    fread(&gid_shell, sizeof(int), 1, in_file);
  }


  fread(&wd->jz_i, sizeof(int), 1, in_file);
  fread(&wd->parity_i, sizeof(int), 1, in_file);
  fread(&wd->w_max_i, sizeof(int), 1, in_file);
  fread(&junk, sizeof(int), 1, in_file);
  fread(&junk, sizeof(int), 1, in_file);
  fread(&junk, sizeof(int), 1, in_file);
  fread(&wd->n_states_i, sizeof(long long int), 1, in_file);
  if (VERBOSE) {printf("The model space has %d shells for a total basis size of %d\n", wd->n_shells, wd->n_states_i);}
  fread(&wd->n_eig_i, sizeof(int), 1, in_file);
  if (VERBOSE) {printf("Initial state file contains %d eigenstates\n", wd->n_eig_i);}
  if (wd->parity_i == '+') {if (VERBOSE) {printf("The initial state basis is positive parity\n");}} 
  else if (wd->parity_i == '-') {if (VERBOSE) {printf("The initial state basis is negative parity\n");}} 
  else if (wd->parity_i == '0') {if (VERBOSE) {printf("The initial state basis contains both positive and negative parity states\n");}} 
  else {printf("Error in initial state parity: %c\n", wd->parity_i); exit(0);}
  if (wd->w_max_i > 0) {
    if (VERBOSE) {
      printf("The max initial state excitation is w_max = %d\n", wd->w_max_i);
    }
  } else if (wd->w_max_i < 0) {
    printf("Error: Max excitation w_max < 0. Please define a truncation scheme where w_max is positive\n"); 
    exit(0);
  }

  wd->e_nuc_i =  (float*) malloc(sizeof(float)*wd->n_eig_i);
  wd->j_nuc_i =  (float*) malloc(sizeof(float)*wd->n_eig_i);
  wd->t_nuc_i =  (float*) malloc(sizeof(float)*wd->n_eig_i);
  unsigned int wave_size = wd->n_states_i*wd->n_eig_i;
  wd->bc_i = (float*) calloc(wave_size, sizeof(float));
  
  if (VERBOSE) {printf("Reading in initial state wavefunction coefficients\n");}
  for (unsigned int i = 0; i < wd->n_eig_i; i++) {
    fread(&junk, sizeof(int), 1, in_file);
    int vec_index;
    double total = 0.0;
    float j_nuc, t_nuc;
    fread(&vec_index, sizeof(int), 1, in_file);
    fread(&wd->e_nuc_i[i], sizeof(float), 1, in_file);
    fread(&j_nuc, sizeof(float), 1, in_file);
    fread(&t_nuc, sizeof(float), 1, in_file);
    t_nuc = -0.5 + 0.5*sqrt(1 + 4.0*t_nuc);
    int ij_nuc = round(2*j_nuc);
    int it_nuc = round(2*t_nuc);
    wd->j_nuc_i[i] = ij_nuc/2.0;
    wd->t_nuc_i[i] = it_nuc/2.0;
    for (unsigned int j = 0; j < wd->n_states_i; j++) {
      unsigned int index = i + wd->n_eig_i*j;
      //fread(&wd->bc_i[i + wd->n_eig_i*j], sizeof(float), 1, in_file);
     // total += pow(wd->bc_i[i + wd->n_eig_i*j], 2);
      fread(&wd->bc_i[index], sizeof(float), 1, in_file);
      total += pow(wd->bc_i[index], 2);

    }
    if (fabs(total -1.0) > pow(10, -6)) {printf("State %d not normalized: norm = %g\n", i, total); exit(0);}
  }
  fclose(in_file);
  if (VERBOSE) {printf("Done.\n");
  printf("Reading in initial state basis\n");}
  in_file = fopen(basis_file_initial, "rb");

  fread(&junk, sizeof(int), 1, in_file);
  fread(&junk, sizeof(int), 1, in_file);
  fread(&junk, sizeof(int), 1, in_file);
  fread(&vec_offset, sizeof(int), 1, in_file);
  fseek(in_file, vec_offset + 16, SEEK_SET);
  for (int i = 0; i < wd->n_shells; i++) {
    fread(&junk, sizeof(int), 1, in_file);
    fread(&wd->n_shell[i], sizeof(int), 1, in_file);
    fread(&wd->l_shell[i], sizeof(int), 1, in_file);
    fread(&wd->j_shell[i], sizeof(int), 1, in_file);
    fread(&wd->jz_shell[i], sizeof(int), 1, in_file);
    fread(&wd->w_shell[i], sizeof(int), 1, in_file);
  }
  for (int i = 0; i < wd->n_shells; i++) {
    int n_shell, j_shell, l_shell, jz_shell, w_shell;
    fread(&junk, sizeof(int), 1, in_file);
    fread(&n_shell, sizeof(int), 1, in_file);
    fread(&l_shell, sizeof(int), 1, in_file);
    fread(&j_shell, sizeof(int), 1, in_file);
    fread(&jz_shell, sizeof(int), 1, in_file);
    fread(&w_shell, sizeof(int), 1, in_file);
  }
  
  wd->wh_hash_i = (wh_list**) calloc(wd->n_sds_p_i*HASH_SIZE, sizeof(wh_list*));
  unsigned int pp0 = n_choose_k(wd->n_shells, wd->n_proton_i);
  unsigned int pn0 = n_choose_k(wd->n_shells, wd->n_neutron_i);
  
  for (int i = 0; i < wd->n_states_i; i++) {
    int in = 0;
    int ip = 0;
    unsigned int pp = pp0;
    unsigned int pn = pn0;
    for (int j = 0; j < wd->n_data; j++) {
      int i_state;
      fread(&i_state, sizeof(int), 1, in_file);
      if (i_state <= wd->n_shells) {
        pp -= n_choose_k(wd->n_shells - i_state, wd->n_proton_i - ip);
        ip++;
      } else {
        pn -= n_choose_k(2*(wd->n_shells) - i_state, wd->n_neutron_i - in);
        in++;
      }
    }
    //printf("%u, %u\n", pp, pn);
    unsigned int p_hash = pp + wd->n_sds_p_i*(pn % HASH_SIZE);
    
    if (wd->wh_hash_i[p_hash] == NULL) {
      wd->wh_hash_i[p_hash] = create_wh_node(pp, pn, i, NULL);
    } else {
      wh_append(wd->wh_hash_i[p_hash], pp, pn, i);
    }
  }
  if (VERBOSE) {printf("Done.\n");}
  fclose(in_file);
  if (VERBOSE) {printf("\n");}
  if (strcmp(wfn_file_initial, wfn_file_final) == 0) {
    if (VERBOSE) {printf("Initial and final bases are identical.\n");}
    wd->same_basis = 1;
    wd->n_proton_f = wd->n_proton_i;
    wd->n_neutron_f = wd->n_neutron_i;
    wd->n_states_f = wd->n_states_i;
    wd->n_sds_p_f = wd->n_sds_p_i;
    wd->n_sds_n_f = wd->n_sds_n_i;
    wd->jz_f = wd->jz_i;
    wd->n_eig_f = wd->n_eig_i;    
    wd->e_nuc_f =  (float*) malloc(sizeof(float)*wd->n_eig_f);
    wd->j_nuc_f =  (float*) malloc(sizeof(float)*wd->n_eig_f);
    wd->t_nuc_f =  (float*) malloc(sizeof(float)*wd->n_eig_f);
    wd->e_nuc_f = wd->e_nuc_i;
    wd->j_nuc_f = wd->j_nuc_i;
    wd->t_nuc_f = wd->t_nuc_i;
    wd->bc_f = wd->bc_i;
    wd->wh_hash_f = wd->wh_hash_i;

  } else {
    if (VERBOSE) {printf("Initial and final bases are different.\n");}
    wd->same_basis = 0;
    in_file = fopen(wfn_file_final, "rb");
    fread(&junk, sizeof(int), 1, in_file);
    fread(&junk, sizeof(int), 1, in_file);
    fread(&vec_offset, sizeof(int), 1, in_file);
    fread(&junk, sizeof(int), 1, in_file);
    fread(&junk, sizeof(int), 1, in_file);
    fread(&junk, sizeof(int), 1, in_file);
    fread(&junk, sizeof(int), 1, in_file);
    fread(&junk, sizeof(int), 1, in_file);
    fread(&wd->n_proton_f, sizeof(int), 1, in_file);
    fread(&wd->n_neutron_f, sizeof(int), 1, in_file);
    if (VERBOSE) {printf("Final state contains %d protons and %d neutrons\n", wd->n_proton_f, wd->n_neutron_f);}
    wd->n_sds_p_f = get_num_sds(wd->n_shells, wd->n_proton_f);
    wd->n_sds_n_f = get_num_sds(wd->n_shells, wd->n_neutron_f);
    if (VERBOSE) {printf("Final proton SDs: %d Final neutron SDs: %d\n", wd->n_sds_p_f, wd->n_sds_n_f);}
    fread(&junk, sizeof(int), 1, in_file);
    fread(&n_proton_orbs, sizeof(int), 1, in_file);
    fread(&n_neutron_orbs, sizeof(int), 1, in_file);
    for (int i = 0; i < n_proton_orbs; i++) {
      int j_orb, pi_orb, w_orb;
      fread(&wd->n_orb[i], sizeof(int), 1, in_file);
      fread(&j_orb, sizeof(int), 1, in_file);
      wd->j_orb[i] = j_orb/2.0;
      fread(&wd->l_orb[i], sizeof(int), 1, in_file);
      fread(&pi_orb, sizeof(int), 1, in_file);
      fread(&w_orb, sizeof(int), 1, in_file);
    }
    for (int i = 0; i < n_neutron_orbs; i++) {
      int n_orb, j_orb, l_orb, pi_orb, w_orb;
      fread(&n_orb, sizeof(int), 1, in_file);
      fread(&j_orb, sizeof(int), 1, in_file);
      fread(&l_orb, sizeof(int), 1, in_file);
      fread(&pi_orb, sizeof(int), 1, in_file);
      fread(&w_orb, sizeof(int), 1, in_file);
    }
    fread(&n_sps_p, sizeof(int), 1, in_file);
    fread(&n_sps_n, sizeof(int), 1, in_file);
    for (int i = 0; i < n_sps_p; i++) {
    int n_shell, j_shell, jz_shell, l_shell, w_shell, pi_shell, id_shell, gid_shell;

      fread(&n_shell, sizeof(int), 1, in_file);
      fread(&j_shell, sizeof(int), 1, in_file);
      fread(&jz_shell, sizeof(int), 1, in_file);
      fread(&l_shell, sizeof(int), 1, in_file);
      fread(&w_shell, sizeof(int), 1, in_file);
      fread(&pi_shell, sizeof(int), 1, in_file);
      fread(&id_shell, sizeof(int), 1, in_file);
      fread(&gid_shell, sizeof(int), 1, in_file);
    }
    for (int i = 0; i < n_sps_n; i++) {
      int n_shell, j_shell, jz_shell, l_shell, w_shell, pi_shell, id_shell, gid_shell;

      fread(&n_shell, sizeof(int), 1, in_file);
      fread(&j_shell, sizeof(int), 1, in_file);
      fread(&jz_shell, sizeof(int), 1, in_file);
      fread(&l_shell, sizeof(int), 1, in_file);
      fread(&w_shell, sizeof(int), 1, in_file);
      fread(&pi_shell, sizeof(int), 1, in_file);
      fread(&id_shell, sizeof(int), 1, in_file);
      fread(&gid_shell, sizeof(int), 1, in_file);
    }

    fread(&wd->jz_f, sizeof(int), 1, in_file);
    fread(&wd->parity_f, sizeof(int), 1, in_file);
    fread(&wd->w_max_f, sizeof(int), 1, in_file);
    fread(&junk, sizeof(int), 1, in_file);
    fread(&junk, sizeof(int), 1, in_file);
    fread(&junk, sizeof(int), 1, in_file);
    fread(&wd->n_states_f, sizeof(long long int), 1, in_file);
    if (VERBOSE) {printf("The model space has %d shells for a total basis size of %d\n", wd->n_shells, wd->n_states_f);}
    fread(&wd->n_eig_f, sizeof(int), 1, in_file);
    if (VERBOSE) {printf("Final state file contains %d eigenstates\n", wd->n_eig_f);}
 
    if (wd->parity_f == '+') {if (VERBOSE) {printf("The final state basis is positive parity\n");}} 
    else if (wd->parity_f == '-') {if (VERBOSE) {printf("The final state basis is negative parity\n");}} 
    else if (wd->parity_f == '0') {if (VERBOSE) {printf("The final state basis contains both positive and negative parity states\n");}} 
    else {printf("Error in final state parity: %c\n", wd->parity_f); exit(0);}
    if (wd->w_max_f > 0) {
      if (VERBOSE) {printf("The max final state excitation is w_max = %d\n", wd->w_max_f);}
    } else if (wd->w_max_f < 0) {
      printf("Error: Max excitation w_max < 0. Please define a truncation scheme where w_max > 0.\n");
      exit(0);
    }

    wd->e_nuc_f = (float*) malloc(sizeof(float)*wd->n_eig_f);
    wd->j_nuc_f = (float*) malloc(sizeof(float)*wd->n_eig_f);
    wd->t_nuc_f = (float*) malloc(sizeof(float)*wd->n_eig_f);
    wave_size = wd->n_states_f*wd->n_eig_f;
    wd->bc_f = (float*) calloc(wave_size, sizeof(float));

    //wd->bc_f = malloc(sizeof(float)*wd->n_states_f*wd->n_eig_f);

    if (VERBOSE) {printf("Reading in final state wavefunction coefficients\n");}
    for (unsigned int i = 0; i < wd->n_eig_f; i++) {
      fread(&junk, sizeof(int), 1, in_file);
      int vec_index;
      float j_nuc, t_nuc;
      double total = 0.0;
      fread(&vec_index, sizeof(int), 1, in_file);
      fread(&wd->e_nuc_f[i], sizeof(float), 1, in_file);
      fread(&j_nuc, sizeof(float), 1, in_file);
      fread(&t_nuc, sizeof(float), 1, in_file);
      t_nuc = -0.5 + 0.5*sqrt(1 + 4.0*t_nuc);
      int ij_nuc = round(2*j_nuc);
      int it_nuc = round(2*t_nuc);
      wd->j_nuc_f[i] = ij_nuc/2.0;
      wd->t_nuc_f[i] = it_nuc/2.0;
      int i_wfn = 0;
      for (unsigned int j = 0; j < wd->n_states_f; j++) {
        fread(&wd->bc_f[i + wd->n_eig_f*j], sizeof(float), 1, in_file);
        total += pow(wd->bc_f[i + wd->n_eig_f*j], 2.0);
	wd->bc_f[i + wd->n_eig_f*j] = 0;
      }
      if (fabs(total - 1.0) > pow(10, -6)) {printf("State %d is not normalized: norm = %g\n", i, total); exit(0);}
    }
    fclose(in_file);
    if (VERBOSE) {printf("Done.\n");
    printf("Reading in final state basis\n");}
    in_file = fopen(basis_file_final, "rb");
    fread(&junk, sizeof(int), 1, in_file);
    fread(&junk, sizeof(int), 1, in_file);
    fread(&junk, sizeof(int), 1, in_file);
    fread(&vec_offset, sizeof(int), 1, in_file);
    fseek(in_file, vec_offset + 16, SEEK_SET);
    for (int i = 0; i < 2*wd->n_shells; i++) {
      int n_shell, j_shell, l_shell, jz_shell, w_shell;
      fread(&junk, sizeof(int), 1, in_file);
      fread(&n_shell, sizeof(int), 1, in_file);
      fread(&l_shell, sizeof(int), 1, in_file);
      fread(&j_shell, sizeof(int), 1, in_file);
      fread(&jz_shell, sizeof(int), 1, in_file);
      fread(&w_shell, sizeof(int), 1, in_file);
      //printf("%d, %d, %d, %d, %d\n", n_shell, l_shell, j_shell, jz_shell, w_shell);
    }
    wd->wh_hash_f = (wh_list**) calloc(wd->n_sds_p_f*HASH_SIZE, sizeof(wh_list*));
    int pp0 = n_choose_k(wd->n_shells, wd->n_proton_f);
    int pn0 = n_choose_k(wd->n_shells, wd->n_neutron_f);
    for (int i = 0; i < wd->n_states_f; i++) {
      int in = 0;
      int ip = 0;
      unsigned int pp = pp0;
      unsigned int pn = pn0;
      for (int j = 0; j < wd->n_data; j++) {
        int i_state;
        fread(&i_state, sizeof(int), 1, in_file);
        if (i_state <= wd->n_shells) {
          pp -= n_choose_k(wd->n_shells - i_state, wd->n_proton_f - ip);
          ip++;
        } else {
          pn -= n_choose_k(2*(wd->n_shells) - i_state, wd->n_neutron_f - in);
          in++;
        }
      }
      unsigned int p_hash = pp + wd->n_sds_p_f*(pn % HASH_SIZE);
      if (wd->wh_hash_f[p_hash] == NULL) {
        wd->wh_hash_f[p_hash] = create_wh_node(pp, pn, i, NULL);
      } else {
        wh_append(wd->wh_hash_f[p_hash], pp, pn, i);
      }
    }
    if (VERBOSE) {printf("Done.\n");}
    fclose(in_file);

  }

  return wd;
}

void modify_binary_wfn_strength(char *wfn_file, float* bc_f) {

  wfnData *wd = malloc(sizeof(*wd));
  FILE *in_file;
  // Read in initial wavefunction data
  if (VERBOSE) {printf("\n");
  printf("Opening initial state wfn file\n");}
  in_file = fopen(wfn_file, "r+b");
  int junk;
  fread(&junk, sizeof(int), 1, in_file);
  fread(&junk, sizeof(int), 1, in_file);
  int vec_offset;
  char buffer[100];
  fread(&vec_offset, sizeof(int), 1, in_file);
  fread(&junk, sizeof(int), 1, in_file);
  fread(&junk, sizeof(int), 1, in_file);
  fread(&junk, sizeof(int), 1, in_file);
  fread(&junk, sizeof(int), 1, in_file);
  fread(&junk, sizeof(int), 1, in_file);
  fread(&wd->n_proton_i, sizeof(int), 1, in_file);
  fread(&wd->n_neutron_i, sizeof(int), 1, in_file);
  wd->n_data = wd->n_proton_i + wd->n_neutron_i;
  if (VERBOSE) {printf("Initial state contains %d protons and %d neutrons\n", wd->n_proton_i, wd->n_neutron_i);}

  fread(&junk, sizeof(int), 1, in_file);
//  fread(&junk, sizeof(int), 1, in_file);
  int n_proton_orbs, n_neutron_orbs;
  fread(&n_proton_orbs, sizeof(int), 1, in_file);
  fread(&n_neutron_orbs, sizeof(int), 1, in_file);
  wd->n_orbits = n_proton_orbs;
  if (VERBOSE) {printf("Model space contains %d orbitals\n", n_neutron_orbs);}
  wd->n_orb = (int*) malloc(sizeof(int)*wd->n_orbits);
  wd->l_orb = (int*) malloc(sizeof(int)*wd->n_orbits);
  wd->w_orb = (int*) malloc(sizeof(int)*wd->n_orbits);
  wd->j_orb = (float*) malloc(sizeof(float)*wd->n_orbits);
  for (int i = 0; i < n_proton_orbs; i++) {
    int j_orb, pi_orb, w_orb;
    fread(&wd->n_orb[i], sizeof(int), 1, in_file);
    fread(&j_orb, sizeof(int), 1, in_file);
    wd->j_orb[i] = j_orb/2.0;
    fread(&wd->l_orb[i], sizeof(int), 1, in_file);
    fread(&pi_orb, sizeof(int), 1, in_file);
    fread(&wd->w_orb[i], sizeof(int), 1, in_file);
  }
  for (int i = 0; i < n_neutron_orbs; i++) {
    int n_orb, j_orb, l_orb, pi_orb, w_orb;
    fread(&n_orb, sizeof(int), 1, in_file);
    fread(&j_orb, sizeof(int), 1, in_file);
    fread(&l_orb, sizeof(int), 1, in_file);
    fread(&pi_orb, sizeof(int), 1, in_file);
    fread(&w_orb, sizeof(int), 1, in_file);
  }
  int n_sps_p, n_sps_n;
  fread(&n_sps_p, sizeof(int), 1, in_file);
  fread(&n_sps_n, sizeof(int), 1, in_file);
  wd->n_shells = n_sps_p;
  if (VERBOSE) {printf("Number of single particle states: %d\n", wd->n_shells);}
  wd->n_shell = (int*) malloc(sizeof(int)*wd->n_shells);
  wd->j_shell = (int*) malloc(sizeof(int)*wd->n_shells);
  wd->l_shell = (int*) malloc(sizeof(int)*wd->n_shells);
  wd->jz_shell = (int*) malloc(sizeof(int)*wd->n_shells);
  wd->w_shell = (int*) malloc(sizeof(int)*wd->n_shells);
  wd->n_sds_p_i = get_num_sds(wd->n_shells, wd->n_proton_i);
  wd->n_sds_n_i = get_num_sds(wd->n_shells, wd->n_neutron_i);
  if (VERBOSE) {printf("Initial proton SDs: %d Initial neutron SDs: %d\n", wd->n_sds_p_i, wd->n_sds_n_i);}
  for (int i = 0; i < n_sps_p; i++) {
    int n_shell, l_shell, j_shell, jz_shell, w_shell, pi_shell, id_shell, gid_shell;
    fread(&n_shell, sizeof(int), 1, in_file);
    fread(&j_shell, sizeof(int), 1, in_file);
    fread(&jz_shell, sizeof(int), 1, in_file);
    fread(&l_shell, sizeof(int), 1, in_file);
    fread(&w_shell, sizeof(int), 1, in_file);
    fread(&pi_shell, sizeof(int), 1, in_file);
    fread(&id_shell, sizeof(int), 1, in_file);
    fread(&gid_shell, sizeof(int), 1, in_file);
  }
  for (int i = 0; i < n_sps_n; i++) {
    int n_shell, j_shell, jz_shell, l_shell, w_shell, pi_shell, id_shell, gid_shell;

    fread(&n_shell, sizeof(int), 1, in_file);
    fread(&j_shell, sizeof(int), 1, in_file);
    fread(&jz_shell, sizeof(int), 1, in_file);
    fread(&l_shell, sizeof(int), 1, in_file);
    fread(&w_shell, sizeof(int), 1, in_file);
    fread(&pi_shell, sizeof(int), 1, in_file);
    fread(&id_shell, sizeof(int), 1, in_file);
    fread(&gid_shell, sizeof(int), 1, in_file);
  }


  fread(&wd->jz_i, sizeof(int), 1, in_file);
  fread(&wd->parity_i, sizeof(int), 1, in_file);
  fread(&wd->w_max_i, sizeof(int), 1, in_file);
  fread(&junk, sizeof(int), 1, in_file);
  fread(&junk, sizeof(int), 1, in_file);
  fread(&junk, sizeof(int), 1, in_file);
  fread(&wd->n_states_i, sizeof(long long int), 1, in_file);
  if (VERBOSE) {printf("The model space has %d shells for a total basis size of %d\n", wd->n_shells, wd->n_states_i);}
  fread(&wd->n_eig_i, sizeof(int), 1, in_file);
  if (VERBOSE) {printf("Initial state file contains %d eigenstates\n", wd->n_eig_i);}
  if (wd->parity_i == '+') {if (VERBOSE) {printf("The initial state basis is positive parity\n");}} 
  else if (wd->parity_i == '-') {if (VERBOSE) {printf("The initial state basis is negative parity\n");}} 
  else if (wd->parity_i == '0') {if (VERBOSE) {printf("The initial state basis contains both positive and negative parity states\n");}} 
  else {printf("Error in initial state parity: %c\n", wd->parity_i); exit(0);}
  if (wd->w_max_i > 0) {
    if (VERBOSE) {
      printf("The max initial state excitation is w_max = %d\n", wd->w_max_i);
    }
  } else if (wd->w_max_i < 0) {
    printf("Error: Max excitation w_max < 0. Please define a truncation scheme where w_max is positive\n"); 
    exit(0);
  }

  wd->e_nuc_i =  (float*) malloc(sizeof(float)*wd->n_eig_i);
  wd->j_nuc_i =  (float*) malloc(sizeof(float)*wd->n_eig_i);
  wd->t_nuc_i =  (float*) malloc(sizeof(float)*wd->n_eig_i);
  unsigned int wave_size = wd->n_states_i*wd->n_eig_i;
  wd->bc_i = (float*) calloc(wave_size, sizeof(float));
  
  if (VERBOSE) {printf("Reading in initial state wavefunction coefficients\n");}
  for (unsigned int i = 0; i < wd->n_eig_i; i++) {
    fread(&junk, sizeof(int), 1, in_file);
    int vec_index;
    double total = 0.0;
    float j_nuc, t_nuc;
    fread(&vec_index, sizeof(int), 1, in_file);
    fread(&wd->e_nuc_i[i], sizeof(float), 1, in_file);
    fread(&j_nuc, sizeof(float), 1, in_file);
    fread(&t_nuc, sizeof(float), 1, in_file);
    t_nuc = -0.5 + 0.5*sqrt(1 + 4.0*t_nuc);
    int ij_nuc = round(2*j_nuc);
    int it_nuc = round(2*t_nuc);
    wd->j_nuc_i[i] = ij_nuc/2.0;
    wd->t_nuc_i[i] = it_nuc/2.0;
    for (unsigned int j = 0; j < wd->n_states_i; j++) {
      unsigned int index = i + wd->n_eig_i*j;
      //fread(&wd->bc_i[i + wd->n_eig_i*j], sizeof(float), 1, in_file);
     // total += pow(wd->bc_i[i + wd->n_eig_i*j], 2);
      fwrite(&bc_f[index], sizeof(float), 1, in_file);
  //    total += pow(wd->bc_i[index], 2);

    }
  //  if (fabs(total -1.0) > pow(10, -6)) {printf("State %d not normalized: norm = %g\n", i, total); exit(0);}
  }
  fclose(in_file);

  return;
}


/* eigen_list* eigen_append(eigen_list* head, int eig_i, int eig_f) {
  if (head->next == NULL) {
    eigen_list* new_node = create_eigen_node(eig_i, eig_f, NULL);
    head->next = new_node;
  } else {
    eigen_list* next = head->next;
    eigen_list* new_node = create_eigen_node(eig_i, eig_f, next);
    head->next = new_node;
  }

  return head;
}
*/
