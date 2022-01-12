#ifndef DENSITY_H
#define DENSITY_H
#include "file_io.h"

void one_body_density(speedParams* sp);
 
void two_body_density(speedParams* sp);



/*   -------------------------
 *   Two-body routines       
 *   -------------------------
 */

// Two-body node trace routines

void trace_two_body_nodes_dmt0_40(int a, int b, int c, int d, int num_mj_p_i, float mj_min_p_i, int num_mj_n_i, float mj_min_n_i, int n_sds_int2, int* p2_array_f, sd_list** p1_list_i, wf_list** n0_list_i, wfnData* wd, int i_op, eigen_list* transition, double* density);

void trace_two_body_nodes_dmt0_22(int a, int b, int c, int d, int num_mj_p_i, float mj_min_p_i, float mj_min_p_f, float mj_max_p_f, int num_mj_n_i, float mj_min_n_i, float mj_min_n_f, float mj_max_n_f, int n_sds_p_int1, int n_sds_n_int1, sd_list** p1_list_i, sd_list** n1_list_i, int* p1_list_f, int* n1_list_f, wfnData* wd, eigen_list* transition, double* density, int *jz_shell); 

void trace_two_body_nodes_dmt0_22_alt(int a, int b, int c, int d, int num_mj_p_i, float mj_min_p_i, float mj_min_p_f, float mj_max_p_f, int num_mj_n_i, float mj_min_n_i, float mj_min_n_f, float mj_max_n_f, sd_list** p1_list_if, sd_list** n1_list_if, wfnData* wd, eigen_list* transition, double* density, int *jz_shell); 

void trace_two_body_nodes_dmtp1_31(int a, int b, int c, int d, int num_mj_p_i, float mj_min_p_i, int num_mj_n_i, float mj_min_n_i, int n_sds_int2, int* p2_array_f, sd_list** p1_list_i, sd_list** n1_list_i, wfnData* wd, int i_op, eigen_list* transition, double* density);

void trace_two_body_nodes_dmtp1_13(int a, int b, int c, int d, int num_mj_p_i, float mj_min_p_i, int num_mj_n_i, float mj_min_n_i, int n_sds_int1, int* p1_array_f, sd_list** p2_list_i, sd_list** n1_list_f, wfnData* wd, int i_op, eigen_list* transition, double* density);

void trace_two_body_nodes_dmtp2(int a, int b, int c, int d, int num_mj_p_i, float mj_min_p_i, int num_mj_n_i, float mj_min_n_i, sd_list** a2_list_i, sd_list** a2_list_f, wfnData* wd, int i_op, eigen_list* transition, double* density);



// Two-body jump routines

void build_two_body_jumps_dmt0(int n_s, int n_p, float mj_min, float mj_max, int num_mj, int n_sds_i, int n_sds_int1, int n_sds_int2, int*a1_array_f, int* a2_array_f, wf_list** a0_list_i, sd_list** a1_list_i, sd_list** a2_list_i, int* jz_shell, int* l_shell);

void build_two_body_jumps_dmtp1(int n_s, int n_proton_i, int num_mj_p_i, float mj_min_p_i, float mj_max_p_i, int n_proton_f, int num_mj_p_f, float mj_min_p_f, float mj_max_p_f, sd_list** p1_list_i, sd_list** p1_list_f, int* p2_array_f, int n_neutron_i, int num_mj_n_i, float mj_min_n_i, float mj_max_n_i, int n_neutron_f, int num_mj_n_f, float mj_min_n_f, float mj_max_n_f, sd_list** n1_list_i, sd_list** n2_list_i, int* n1_array_f, int* jz_shell, int* l_shell);

void build_two_body_jumps_dmtp2(int n_s, int n_proton, int num_mj_p_i, float mj_min_p_i, float mj_max_p_i, int num_mj_p_f, float mj_min_p_f, float mj_max_p_f, int n_sds_p_f, sd_list** p2_list_f, int n_neutron, int num_mj_n_i, float mj_min_n_i, float mj_max_n_i, int num_mj_n_f, float mj_min_n_f, float mj_max_n_f, int n_sds_n_i, sd_list** n2_list_i, int* jz_shell, int* l_shell);

void build_two_body_jumps_dmt0_alt(int n_s, int n_p, float mj_min, float mj_max, int num_mj, int n_sds_i, int n_sds_int1, int n_sds_int2, int*a1_array_f, int* a2_array_f, wf_list** a0_list_i, sd_list** a1_list_i, sd_list** a2_list_i, sd_list** a2_list_if, int* jz_shell, int* l_shell);

// Two-body jumps routines w/ truncation

void build_two_body_jumps_dmt0_trunc(int n_s, int n_p, float mj_min, float mj_max, int num_mj, int n_sds_i, int n_sds_int1, int n_sds_int2, int*a1_array_f, int* a2_array_f, wf_list** a0_list_i, sd_list** a1_list_i, sd_list** a2_list_i, int* jz_shell, int* l_shell, int* w_shell, int w_max);

void build_two_body_jumps_dmtp1_trunc(int n_s, int n_proton_i, int num_mj_p_i, float mj_min_p_i, float mj_max_p_i, int n_proton_f, int num_mj_p_f, float mj_min_p_f, float mj_max_p_f, sd_list** p1_list_i, sd_list** p1_list_f, int* p2_array_f, int n_neutron_i, int num_mj_n_i, float mj_min_n_i, float mj_max_n_i, int n_neutron_f, int num_mj_n_f, float mj_min_n_f, float mj_max_n_f, sd_list** n1_list_i, sd_list** n2_list_i, int* n1_array_f, int* jz_shell, int* l_shell, int* w_shell, int w_max_p_i, int w_max_p_f, int w_max_n_i, int w_max_n_f);

void build_two_body_jumps_dmtp2_trunc(int n_s, int n_proton, int num_mj_p_i, float mj_min_p_i, float mj_max_p_i, int num_mj_p_f, float mj_min_p_f, float mj_max_p_f, int n_sds_p_f, sd_list** p2_list_f, int n_neutron, int num_mj_n_i, float mj_min_n_i, float mj_max_n_i, int num_mj_n_f, float mj_min_n_f, float mj_max_n_f, int n_sds_n_i, sd_list** n2_list_i, int* jz_shell, int* l_shell, int* w_shell, int w_max_p_i, int w_max_p_f, int w_max_n_i, int w_max_n_f); 

/*   -------------------------
 *   One-body routines       
 *   -------------------------
 */

// One-body node trace routines


void trace_one_body_nodes_dmt0(int a, int b, int num_mj_p_i, float min_mj_p_i, int num_mj_n_i, float min_mj_n_i, int n_sds_int, int* a1_array_f, sd_list** a1_list_i, wf_list** a0_list_i, wfnData* wd, int i_op, eigen_list* transition, double *density);

void trace_one_body_nodes_dmtp1(int a, int b, int num_mj_p_i, float mj_min_p_i, int num_mj_n_i, float mj_min_n_i, sd_list** a1_list_i, sd_list** a1_list_f, wfnData* wd, int i_op, eigen_list *transition, double* density);


// One-body jump routines

void build_one_body_jumps_dmt0(int n_s, int n_p, float mj_min, float mj_max, int num_mj, int n_sds_f, int n_sds_int, int* a1_array_f, wf_list** a0_list_i, sd_list** a1_list_i, int* jz_shell, int* l_shell);

void build_one_body_jumps_dmtp1(int n_s, int n_proton_f, int num_mj_p_i, float mj_min_p_i, float mj_max_p_i, int num_mj_p_f, float mj_min_p_f, float mj_max_p_f, sd_list** p1_list_f, int n_neutron_i, int num_mj_n_i, float mj_min_n_i, float mj_max_n_i, sd_list** n1_list_i, int* jz_shell, int* l_shell);


// One-body jump routines w/ truncation

void build_one_body_jumps_dmt0_trunc(int n_s, int n_p, float mj_min, float mj_max, int num_mj, int n_sds_f, int n_sds_int, int* a1_array_f, wf_list** a0_list_i, sd_list** a1_list_i, int* jz_shell, int* l_shell, int* w_shell, int w_max);

void build_one_body_jumps_dmtp1_trunc(int n_s, int n_proton_f, int num_mj_p_i, float mj_min_p_i, float mj_max_p_i, int num_mj_p_f, float mj_min_p_f, float mj_max_p_f, sd_list** p1_list_f, int n_neutron_i, int num_mj_n_i, float mj_min_n_i, float mj_max_n_i, sd_list** n1_list_i, int* jz_shell, int* l_shell, int* w_shell, int w_max_p_i, int w_max_p_f, int w_max_n_i, int w_max_n_f);

int test_suite();
double compute_2body_J0_T0_sum_rule(char *input_file);
double compute_2body_J0_T1_sum_rule(char *input_file);
double compute_2body_J0_T2_sum_rule(char *input_file);


#endif
