#ifndef FILE_IO_H
#define FILE_IO_H
#include "slater.h"

typedef struct eigen_list
{
  int eig_i;
  int eig_f;
  struct eigen_list *next;
} eigen_list;


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

typedef struct sde_list
{
  int pi, pn;
  int phase;
  int n_quanta;
  struct sde_list *next;
} sde_list;


typedef struct wf_list
{
  unsigned int p;
  struct wf_list *next;
} wf_list;

typedef struct wfe_list
{
  unsigned int p;
  int n_quanta;
  struct wfe_list *next;
} wfe_list;


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
  int parity_i, wmax_i, parity_f, wmax_f;
  unsigned int n_sds_p_i, n_sds_p_f, n_sds_n_i, n_sds_n_f;
  wh_list **wh_hash_i, **wh_hash_f;
  float *bc_i, *bc_f;
  int *n_shell, *l_shell, *j_shell, *jz_shell, *tz_shell, *w_shell;
  int *n_orb, *l_orb, *w_orb;
  float *j_orb;
  float jz_i, jz_f;
  float *e_nuc_i, *e_nuc_f, *j_nuc_i, *j_nuc_f, *t_nuc_i, *t_nuc_f;
  int same_basis;
} wfnData;

typedef struct speedParams
{
  char *initial_file_base, *final_file_base, *out_file_base;
  int n_body;
  int n_trans;
  int j_op, t_op;
  eigen_list *transition_list;
  
} speedParams;


sd_list* create_sd_node(unsigned int pi, unsigned int pf, int phase, sd_list* next);
sd_list* sd_append(sd_list* head, unsigned int pi, unsigned int pf, int phase);
sde_list* create_sde_node(unsigned int pi, unsigned int pf, int phase, int n_quanta, sde_list* next);
sde_list* sde_append(sde_list* head, unsigned int pi, unsigned int pf, int phase, int n_quanta);
wf_list* create_wf_node(unsigned int b, wf_list* next);
wf_list* wf_append(wf_list* head, unsigned int p);
wfe_list* create_wfe_node(unsigned int b, int n_quanta, wfe_list* next);
wfe_list* wfe_append(wfe_list* head, unsigned int p, int n_quanta);
wh_list* create_wh_node(unsigned int pp, unsigned int pn, unsigned int index, wh_list* next);
wh_list* wh_append(wh_list* head, unsigned int pp, unsigned int pn, unsigned int index);

eigen_list* create_eigen_node(int eig_i, int eig_n, eigen_list* next);
eigen_list* eigen_append(eigen_list* head, int eig_i, int eig_f);
wfnData* read_binary_wfn_data(char *wfn_file_initial, char *wfn_file_final, char* basis_file_initial, char *basis_file_final);

speedParams* read_parameter_file(char* parameter_file);
#endif
