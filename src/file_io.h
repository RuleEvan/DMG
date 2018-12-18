#ifndef FILE_IO_H
#define FILE_IO_H
#include "slater.h"

typedef struct sBasisCoeff 
{
  double *wave;
} BasisCoeff;

typedef struct sd_list
{
  int pi, pn;
  int phase;
  struct sd_list *next;
} sd_list;

typedef struct wf_list
{
  unsigned int p;
  struct wf_list *next;
} wf_list;

typedef struct wh_list
{
  unsigned int pn;
  unsigned int pp;
  unsigned int index;
  struct wh_list *next;
} wh_list;

typedef struct wfnData
{
  int n_proton_i, n_proton_f, n_neutron_i, n_neutron_f;
  long long int n_states_i, n_states_f;
  int n_shells, n_orbits, n_data;
  int n_eig_i, n_eig_f;
  unsigned int n_sds_p_i, n_sds_p_f, n_sds_n_i, n_sds_n_f;
  wh_list **wh_hash_i, **wh_hash_f;
  float *bc_i, *bc_f;
  int *n_shell, *l_shell, *j_shell, *jz_shell, *tz_shell;
  int *n_orb, *l_orb;
  float *j_orb;
  float jz_i, jz_f;
  float *e_nuc_i, *e_nuc_f, *j_nuc_i, *j_nuc_f, *t_nuc_i, *t_nuc_f;
  int same_basis;
} wfnData;

sd_list* create_sd_node(int pi, int pf, int phase, sd_list* next);
sd_list* sd_append(sd_list* head, int pi, int pf, int phase);
wf_list* create_wf_node(unsigned int b, wf_list* next);
wf_list* wf_append(wf_list* head, unsigned int p);
wh_list* create_wh_node(unsigned int pp, unsigned int pn, unsigned int index, wh_list* next);
wh_list* wh_append(wh_list* head, unsigned int pp, unsigned int pn, unsigned int index);
wfnData* read_wfn_data(char *wfn_file_initial, char *wfn_file_final, char *orbit_file);
wfnData* read_binary_wfn_data(char *wfn_file_initial, char *wfn_file_final, char* basis_file_initial, char *basis_file_final, char *orbit_file);
#endif
