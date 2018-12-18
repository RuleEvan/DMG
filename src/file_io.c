#include "file_io.h"

wfnData* read_binary_wfn_data(char *wfn_file_initial, char *wfn_file_final, char *basis_file_initial, char *basis_file_final, char *orbit_file) {
  wfnData *wd = malloc(sizeof(*wd));
  FILE *in_file;
  // Read in initial wavefunction data
  printf("Opening file\n");
  in_file = fopen(wfn_file_initial, "rb");
  int junk;
  fread(&junk, sizeof(int), 1, in_file);
  fread(&junk, sizeof(int), 1, in_file);
  int vec_offset;
  char buffer[100];
  fread(&vec_offset, sizeof(int), 1, in_file);
  fread(&junk, sizeof(int), 5, in_file);
  fread(&wd->n_proton_i, sizeof(int), 1, in_file);
  fread(&wd->n_neutron_i, sizeof(int), 1, in_file);
  wd->n_data = wd->n_proton_i + wd->n_neutron_i;
  printf("Initial state contains %d protons and %d neutrons\n", wd->n_proton_i, wd->n_neutron_i);

  fread(&junk, sizeof(int), 1, in_file);
//  fread(&junk, sizeof(int), 1, in_file);
  int n_proton_orbs, n_neutron_orbs;
  fread(&n_proton_orbs, sizeof(int), 1, in_file);
  fread(&n_neutron_orbs, sizeof(int), 1, in_file);
  wd->n_orbits = n_proton_orbs;
  printf("Model space contains %d orbitals\n", n_neutron_orbs);
  wd->n_orb = (int*) malloc(sizeof(int)*wd->n_orbits);
  wd->l_orb = (int*) malloc(sizeof(int)*wd->n_orbits);
  wd->j_orb = (float*) malloc(sizeof(float)*wd->n_orbits);
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
  int n_sps_p, n_sps_n;
  fread(&n_sps_p, sizeof(int), 1, in_file);
  fread(&n_sps_n, sizeof(int), 1, in_file);
  wd->n_shells = n_sps_p;
  wd->n_shell = (int*) malloc(sizeof(int)*wd->n_shells);
  wd->j_shell = (int*) malloc(sizeof(int)*wd->n_shells);
  wd->l_shell = (int*) malloc(sizeof(int)*wd->n_shells);
  wd->jz_shell = (int*) malloc(sizeof(int)*wd->n_shells);
  wd->n_sds_p_i = get_num_sds(wd->n_shells, wd->n_proton_i);
  wd->n_sds_n_i = get_num_sds(wd->n_shells, wd->n_neutron_i);
  printf("Initial proton SDs: %d Initial neutron SDs: %d\n", wd->n_sds_p_i, wd->n_sds_n_i);
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
  fread(&junk, sizeof(int), 1, in_file);
  fread(&junk, sizeof(int), 1, in_file);
  fread(&junk, sizeof(int), 1, in_file);
  fread(&junk, sizeof(int), 1, in_file);
  fread(&junk, sizeof(int), 1, in_file);
  fread(&wd->n_states_i, sizeof(long long int), 1, in_file);
  printf("The model space has %d shells for a total basis size of %d\n", wd->n_shells, wd->n_states_i);
  fread(&wd->n_eig_i, sizeof(int), 1, in_file);
  printf("Initial state file contains %d eigenstates\n", wd->n_eig_i);
  
  wd->e_nuc_i =  (float*) malloc(sizeof(float)*wd->n_eig_i);
  wd->j_nuc_i =  (float*) malloc(sizeof(float)*wd->n_eig_i);
  wd->t_nuc_i =  (float*) malloc(sizeof(float)*wd->n_eig_i);
  wd->bc_i = malloc(sizeof(float)*wd->n_states_i*wd->n_eig_i);

  printf("Reading in initial state wavefunction coefficients\n");
  for (int i = 0; i < wd->n_eig_i; i++) {
    fread(&junk, sizeof(int), 1, in_file);
    int vec_index;
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
    for (int j = 0; j < wd->n_states_i; j++) {
      fread(&wd->bc_i[i + wd->n_eig_i*j], sizeof(float), 1, in_file);
    }
  }
  fclose(in_file);
  printf("Done.\n");
  printf("Reading in initial state basis\n");
  in_file = fopen(basis_file_initial, "rb");
  fread(&junk, sizeof(int), 1, in_file);
  fread(&junk, sizeof(int), 1, in_file);
  fread(&junk, sizeof(int), 1, in_file);
  fread(&vec_offset, sizeof(int), 1, in_file);
  fseek(in_file, vec_offset + 16, SEEK_SET);
  for (int i = 0; i < wd->n_shells; i++) {
    int w_shell;
    fread(&junk, sizeof(int), 1, in_file);
    fread(&wd->n_shell[i], sizeof(int), 1, in_file);
    fread(&wd->l_shell[i], sizeof(int), 1, in_file);
    fread(&wd->j_shell[i], sizeof(int), 1, in_file);
    fread(&wd->jz_shell[i], sizeof(int), 1, in_file);
    fread(&w_shell, sizeof(int), 1, in_file);
//    printf("%d, %d, %d, %d, %d\n", i+1, wd->n_shell[i], wd->l_shell[i], wd->j_shell[i], wd->jz_shell[i]);
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
  int *p_orbitals = (int*) malloc(sizeof(int)*wd->n_proton_i);
  int *n_orbitals = (int*) malloc(sizeof(int)*wd->n_neutron_i);
  for (int i = 0; i < wd->n_states_i; i++) {
    int in = 0;
    int ip = 0;
    for (int j = 0; j < wd->n_data; j++) {
      int i_state;
      fread(&i_state, sizeof(int), 1, in_file);
      if (i_state <= wd->n_shells) {
        p_orbitals[ip] = i_state;
        ip++;
      } else {
        n_orbitals[in] = i_state - wd->n_shells;
        in++;
      }
    }
    if ((p_orbitals[0] == 1) && (p_orbitals[1] == 10) && (p_orbitals[2] == 18) && (p_orbitals[3] == 21) && (n_orbitals[0] == 7) && (n_orbitals[1] == 9) && (n_orbitals[2] == 18) && (n_orbitals[3] == 20)) {exit(0);} 
    unsigned int pp = p_step(wd->n_shells, wd->n_proton_i, p_orbitals);
    unsigned int pn = p_step(wd->n_shells, wd->n_neutron_i, n_orbitals);
    unsigned int p_hash = pp + wd->n_sds_p_i*(pn % HASH_SIZE);
    if (wd->wh_hash_i[p_hash] == NULL) {
      wd->wh_hash_i[p_hash] = create_wh_node(pp, pn, i, NULL);
    } else {
      wh_append(wd->wh_hash_i[p_hash], pp, pn, i);
    }
  }
  free(p_orbitals);
  free(n_orbitals);
  printf("Done.\n");
  fclose(in_file);
  
  if (strcmp(wfn_file_initial, wfn_file_final) == 0) {
    printf("Initial and final states are identical\n");
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
    printf("Initial and final states differ\n");
    wd->same_basis = 0;
    in_file = fopen(wfn_file_final, "rb");
    fread(&junk, sizeof(int), 1, in_file);
    fread(&junk, sizeof(int), 1, in_file);
    fread(&vec_offset, sizeof(int), 1, in_file);
    fread(&junk, sizeof(int), 5, in_file);
    fread(&wd->n_proton_f, sizeof(int), 1, in_file);
    fread(&wd->n_neutron_f, sizeof(int), 1, in_file);
    printf("Final state contains %d protons and %d neutrons\n", wd->n_proton_f, wd->n_neutron_f);
    wd->n_sds_p_f = get_num_sds(wd->n_shells, wd->n_proton_f);
    wd->n_sds_n_f = get_num_sds(wd->n_shells, wd->n_neutron_f);
    printf("Final proton SDs: %d Final neutron SDs: %d\n", wd->n_sds_p_f, wd->n_sds_n_f);
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
    fread(&junk, sizeof(int), 1, in_file);
    fread(&junk, sizeof(int), 1, in_file);
    fread(&junk, sizeof(int), 1, in_file);
    fread(&junk, sizeof(int), 1, in_file);
    fread(&junk, sizeof(int), 1, in_file);
    fread(&wd->n_states_f, sizeof(long long int), 1, in_file);
    printf("The model space has %d shells for a total basis size of %d\n", wd->n_shells, wd->n_states_f);
    fread(&wd->n_eig_f, sizeof(int), 1, in_file);
    printf("Final state file contains %d eigenstates\n", wd->n_eig_f);
  
    wd->e_nuc_f = (float*) malloc(sizeof(float)*wd->n_eig_f);
    wd->j_nuc_f = (float*) malloc(sizeof(float)*wd->n_eig_f);
    wd->t_nuc_f = (float*) malloc(sizeof(float)*wd->n_eig_f);
    wd->bc_f = malloc(sizeof(float)*wd->n_states_f*wd->n_eig_f);

    printf("Reading in final state wavefunction coefficients\n");
    for (int i = 0; i < wd->n_eig_f; i++) {
      fread(&junk, sizeof(int), 1, in_file);
      int vec_index;
      float j_nuc, t_nuc;
      fread(&vec_index, sizeof(int), 1, in_file);
      fread(&wd->e_nuc_f[i], sizeof(float), 1, in_file);
      fread(&j_nuc, sizeof(float), 1, in_file);
      fread(&t_nuc, sizeof(float), 1, in_file);
      t_nuc = -0.5 + 0.5*sqrt(1 + 4.0*t_nuc);
      int ij_nuc = round(2*j_nuc);
      int it_nuc = round(2*t_nuc);
      wd->j_nuc_f[i] = ij_nuc/2.0;
      wd->t_nuc_f[i] = it_nuc/2.0;
      for (int j = 0; j < wd->n_states_f; j++) {
        fread(&wd->bc_f[i + wd->n_eig_f*j], sizeof(float), 1, in_file);
      }
    }
    fclose(in_file);
    printf("Done.\n");
    printf("Reading in final state basis\n");
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
    }
  
    wd->wh_hash_f = (wh_list**) calloc(wd->n_sds_p_f*HASH_SIZE, sizeof(wh_list*));
    p_orbitals = (int*) malloc(sizeof(int)*wd->n_proton_f);
    n_orbitals = (int*) malloc(sizeof(int)*wd->n_neutron_f);

    for (int i = 0; i < wd->n_states_f; i++) {
      int in = 0;
      int ip = 0;
      for (int j = 0; j < wd->n_proton_f + wd->n_neutron_f; j++) {
        int i_state;
        fread(&i_state, sizeof(int), 1, in_file);
        if (i_state <= wd->n_shells) {
          p_orbitals[ip] = i_state;
          ip++;
        } else {
          n_orbitals[in] = i_state - wd->n_shells;
          in++;
        }
      }
      unsigned int pp = p_step(wd->n_shells, wd->n_proton_f, p_orbitals);
      unsigned int pn;
      if (in == 0) {pn = 1;} else {pn = p_step(wd->n_shells, wd->n_neutron_f, n_orbitals);}
      unsigned int p_hash = pp + wd->n_sds_p_f*(pn % HASH_SIZE);
      if (wd->wh_hash_f[p_hash] == NULL) {
        wd->wh_hash_f[p_hash] = create_wh_node(pp, pn, i, NULL);
      } else {
        wh_append(wd->wh_hash_f[p_hash], pp, pn, i);
      }
    }
    printf("Done.\n");
    free(p_orbitals);
    free(n_orbitals);
    fclose(in_file);

  }

  return wd;
}


wfnData* read_wfn_data(char *wfn_file_initial, char *wfn_file_final, char *orbit_file) {
  wfnData *wd = malloc(sizeof(*wd));
  FILE* in_file;
  // Read in initial wavefunction data
  in_file = fopen(wfn_file_initial, "r");
  fscanf(in_file, "%d\n", &wd->n_proton_i);
  fscanf(in_file, "%d\n", &wd->n_neutron_i);
  wd->n_data = wd->n_proton_i + wd->n_neutron_i;
  char buffer[100];
  fgets(buffer, 100, in_file);
  fgets(buffer, 100, in_file);
  fgets(buffer, 100, in_file);
  fscanf(in_file, "%d", &wd->n_shells);
  wd->n_shells /= 2; // Mod out isospin
  fgets(buffer, 100, in_file);
  fgets(buffer, 100, in_file);
  fscanf(in_file, "%d", &wd->n_states_i);
  printf("Initial state contains %d protons and %d neutrons\n", wd->n_proton_i, wd->n_neutron_i);
  printf("The model space has %d shells for a total basis size of %d\n", wd->n_shells, wd->n_states_i);

  fgets(buffer, 100, in_file);
  fgets(buffer, 100, in_file);
  fscanf(in_file, "%f", &wd->jz_i);
  fgets(buffer, 100, in_file);
  fscanf(in_file, "%d", &wd->n_eig_i);
  printf("Initial state file contains %d eigenstates\n", wd->n_eig_i);
  fgets(buffer, 100, in_file);
  wd->e_nuc_i =  (float*) malloc(sizeof(float)*wd->n_eig_i);
  wd->j_nuc_i =  (float*) malloc(sizeof(float)*wd->n_eig_i);
  wd->t_nuc_i =  (float*) malloc(sizeof(float)*wd->n_eig_i);
  for (int i = 0; i < wd->n_eig_i; i++) {
    fscanf(in_file, "%f", &wd->e_nuc_i[i]);
    fscanf(in_file, "%f", &wd->j_nuc_i[i]);
    fscanf(in_file, "%f", &wd->t_nuc_i[i]);
    wd->j_nuc_i[i] = round(fabs(wd->j_nuc_i[i]));
    wd->t_nuc_i[i] = round(fabs(wd->t_nuc_i[i]));
    if (wd->n_data % 2) {
      wd->j_nuc_i[i] -= 0.5;
      wd->t_nuc_i[i] -= 0.5;
    }

    fgets(buffer, 100, in_file);
  }
  wd->n_shell = (int*) malloc(sizeof(int)*wd->n_shells);
  wd->j_shell = (int*) malloc(sizeof(int)*wd->n_shells);
  wd->l_shell = (int*) malloc(sizeof(int)*wd->n_shells);
  wd->jz_shell = (int*) malloc(sizeof(int)*wd->n_shells);
  //wd->tz_shell = (int*) malloc(sizeof(int)*wd->n_shells);
  wd->n_sds_p_i = get_num_sds(wd->n_shells, wd->n_proton_i);
  wd->n_sds_n_i = get_num_sds(wd->n_shells, wd->n_neutron_i);

  // Read in shell quantum numbers
  for (int i = 0; i < wd->n_shells; i++) {
    fscanf(in_file, "%*d %d %d %d %d %*d\n", &wd->n_shell[i], &wd->l_shell[i], &wd->j_shell[i], &wd->jz_shell[i]);
  }
  // Skip shells with opposite isospin
  for (int i = 0; i < wd->n_shells; i++) {
    fscanf(in_file, "%*d %*d %*d %*d %*d %*d\n");
  }

  printf("Reading in initial state wavefunction coefficients\n");
  wd->wh_hash_i = (wh_list**) calloc(HASH_SIZE, sizeof(wh_list*));
  wd->bc_i = malloc(sizeof(float)*wd->n_states_i*wd->n_eig_i);
  int *p_orbitals = (int*) malloc(sizeof(int)*wd->n_proton_i);
  int *n_orbitals = (int*) malloc(sizeof(int)*wd->n_neutron_i);
  for (int i = 0; i < wd->n_states_i; i++) {
    int in = 0;
    int ip = 0;
    for (int j = 0; j < wd->n_data; j++) {
      int i_orb;
      fscanf(in_file, "%d", &i_orb);
      if (i_orb <= wd->n_shells) {
        p_orbitals[ip] = i_orb;
        ip++;
      } else {
        n_orbitals[in] = i_orb - wd->n_shells;
        in++;
      }
    }
    unsigned int pp = p_step(wd->n_shells, wd->n_proton_i, p_orbitals);
    unsigned int pn = p_step(wd->n_shells, wd->n_neutron_i, n_orbitals);
    unsigned int p_hash = pp + wd->n_sds_p_i*(pn % HASH_SIZE);
    if (wd->wh_hash_i[p_hash] == NULL) {
      wd->wh_hash_i[p_hash] = create_wh_node(pp, pn, i, NULL);
    } else {
      wh_append(wd->wh_hash_i[p_hash], pp, pn, i);
    }

    fgets(buffer, 100, in_file);
    for (int j = 0; j < wd->n_eig_i; j++) {
      fscanf(in_file, "%f\n", &wd->bc_i[j + wd->n_eig_i*i]);
    }
  }
  free(p_orbitals);
  free(n_orbitals);
  fclose(in_file);
  
    

  if (strcmp(wfn_file_initial, wfn_file_final) == 0) {
    printf("Initial and final states are identical\n");
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
    wd->same_basis = 0;
    in_file = fopen(wfn_file_final, "r");
    int n_shells_test;
    fscanf(in_file, "%d\n", &wd->n_proton_f);
    fscanf(in_file, "%d\n", &wd->n_neutron_f);
    wd->n_sds_p_f = get_num_sds(wd->n_shells, wd->n_proton_f);
    wd->n_sds_n_f = get_num_sds(wd->n_shells, wd->n_neutron_f);
    fgets(buffer, 100, in_file);
    fgets(buffer, 100, in_file);
    fgets(buffer, 100, in_file);
    fscanf(in_file, "%d", &n_shells_test);
    fgets(buffer, 100, in_file);
    fgets(buffer, 100, in_file);
    fscanf(in_file, "%d", &wd->n_states_f);
    printf("Final state contains %d protons and %d neutrons\n", wd->n_proton_f, wd->n_neutron_f);
    if (wd->n_shells != n_shells_test/2) {printf("Error: number of shells does not agree between initial and final state model spaces\n"); exit(0);}
    printf("The model space has %d shells for a total basis size of %d\n", wd->n_shells, wd->n_states_f);
    if (wd->n_proton_f + wd->n_neutron_f != wd->n_proton_i + wd->n_neutron_i) {printf("Error: total number of nucleons is not constant\n"); exit(0);}
    fgets(buffer, 100, in_file);
    fgets(buffer, 100, in_file);
    fscanf(in_file, "%f", &wd->jz_f);
    fgets(buffer, 100, in_file);
    fscanf(in_file, "%d", &wd->n_eig_f);
    printf("Final state file contains %d eigenvalues\n", wd->n_eig_f);
    fgets(buffer, 100, in_file);
    wd->e_nuc_f =  (float*) malloc(sizeof(float)*wd->n_eig_f);
    wd->j_nuc_f =  (float*) malloc(sizeof(float)*wd->n_eig_f);
    wd->t_nuc_f =  (float*) malloc(sizeof(float)*wd->n_eig_f);

    for (int i = 0; i < wd->n_eig_f; i++) {
      fscanf(in_file, "%f", &wd->e_nuc_f[i]);
      fscanf(in_file, "%f", &wd->j_nuc_f[i]);
      fscanf(in_file, "%f", &wd->t_nuc_f[i]);
      wd->j_nuc_f[i] = round(fabs(wd->j_nuc_f[i]));
      wd->t_nuc_f[i] = round(fabs(wd->t_nuc_f[i]));
      if (wd->n_data % 2) {
        wd->j_nuc_f[i] -= 0.5;
        wd->t_nuc_f[i] -= 0.5;
      }
      fgets(buffer, 100, in_file);
    }
    for (int i = 0; i < 2*wd->n_shells; i++) {
      fscanf(in_file, "%*d %*d %*d %*d %*d %*d\n");
    }
    wd->wh_hash_f = (wh_list**) calloc(HASH_SIZE, sizeof(wh_list*));
    p_orbitals = (int*) malloc(sizeof(int)*wd->n_proton_f);
    n_orbitals = (int*) malloc(sizeof(int)*wd->n_neutron_f);
    wd->bc_f = malloc(sizeof(float)*wd->n_states_f*wd->n_eig_f);
    for (int i = 0; i < wd->n_states_f; i++) {
      int in = 0;
      int ip = 0;
      for (int j = 0; j < wd->n_data; j++) {
        int i_orb;
        fscanf(in_file, "%d", &i_orb);
        if (i_orb <= wd->n_shells) {
          p_orbitals[ip] = i_orb;
          ip++;
        } else {
          n_orbitals[in] = i_orb - wd->n_shells;
          in++;
        }
      }
      unsigned int pp = p_step(wd->n_shells, wd->n_proton_f, p_orbitals);
      unsigned int pn;
      if (in == 0) {pn = 1;
      } else {
        pn = p_step(wd->n_shells, wd->n_neutron_f, n_orbitals);
      }
      unsigned int p_hash = pp + wd->n_sds_p_f*(pn % HASH_SIZE);
      if (wd->wh_hash_f[p_hash] == NULL) {
        wd->wh_hash_f[p_hash] = create_wh_node(pp, pn, i, NULL);
      } else {
        wh_append(wd->wh_hash_f[p_hash], pp, pn, i);
      }
      fgets(buffer, 100, in_file);
      for (int j = 0; j < wd->n_eig_f; j++) {
        fscanf(in_file, "%f\n", &wd->bc_f[j + wd->n_eig_f*i]);
      }
    }
    free(p_orbitals);
    free(n_orbitals);

    fclose(in_file);
  }

  in_file = fopen(orbit_file, "r");
  fgets(buffer, 100, in_file);
  fscanf(in_file, "%d\n", &wd->n_orbits);
  printf("The model space contains %d orbitals\n", wd->n_orbits);
  wd->n_orb = (int*) malloc(sizeof(int)*wd->n_orbits);
  wd->l_orb = (int*) malloc(sizeof(int)*wd->n_orbits);
  wd->j_orb = (float*) malloc(sizeof(float)*wd->n_orbits);
  for (int i = 0; i < wd->n_orbits; i++) {
    float n_orb_f, l_orb_f;
    fscanf(in_file, "%f %f %f %*d", &n_orb_f, &l_orb_f, &wd->j_orb[i]);
    wd->n_orb[i] = (int) n_orb_f;
    wd->l_orb[i] = (int) l_orb_f;
  }
  fclose(in_file);
  return wd;
}


sd_list* create_sd_node(int pi, int pn, int phase, sd_list* next) {
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

sd_list* sd_append(sd_list* head, int pi, int pn, int phase) {
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



