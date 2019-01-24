#ifndef DENSITY_H
#define DENSITY_H
#include "file_io.h"

void one_body_density(int j_op, int t_op);

void two_body_density(int j_op, int t_op);

double trace_a4_nodes(int a, int b, int c, int d, int num_mj, int n_sds_int2, int* p2_array_f, sd_list** p1_list_i, wf_list** n0_list_i, wfnData* wd, int psi_i, int psi_f, int i_op);

double trace_a22_nodes(int a, int b, int c, int d, int num_mj, sd_list** a2_list_i, sd_list** a2_list_f, wfnData* wd, int psi_i, int psi_f, int i_op);

double trace_a20_nodes(int a, int b, int c, int d, int num_mj, int n_sds_p_int1, int n_sds_n_int1, sd_list** p1_list_i, sd_list** n1_list_i, int* p1_list_f, int* n1_list_f, wfnData* wd, int psi_i, int psi_f); 

void build_two_body_jumps_f(int n_s, int n_p, float mj_min, float mj_max, int num_mj, int n_sds_f, int n_sds_int1, int n_sds_int2, int* a1_array_f, int* a2_array_f, sd_list** a2_list_f, int* jz_shell, int* l_shell);

void build_two_body_jumps_i(int n_s, int n_p, float mj_min, float mj_max, int num_mj, int n_sds_i, wf_list** a0_list_i, sd_list** a1_list_i, sd_list** a2_list_i, int* jz_shell, int* l_shell);

void build_two_body_jumps_i_and_f(int n_s, int n_p, float mj_min, float mj_max, int num_mj, int n_sds_f, int n_sds_int1, int n_sds_int2, int* a1_array_f, int* a2_array_f, wf_list** a0_list_i, sd_list** a1_list_i, sd_list** a2_list_i, int* jz_shell, int* l_shell);

double trace_1body_t0_nodes(int a, int b, int num_mj, int n_sds_int, int* a1_array_f, sd_list** a1_list_i, wf_list** a0_list_i, wfnData* wd, int psi_i, int psi_f, int i_op);

double trace_1body_t2_nodes(int a, int b, int num_mj, sd_list** a1_list_i, sd_list** a1_list_f, wfnData* wd, int psi_i, int psi_f, int i_op);

void build_one_body_jumps_f(int n_s, int n_p, float mj_min, float mj_max, int num_mj, int n_sds_f, int n_sds_int, int* a1_array_f, sd_list** a1_list_f, int* jz_shell, int* l_shell);

void build_one_body_jumps_i(int n_s, int n_p, float mj_min, float mj_max, int num_mj, int n_sds_i, wf_list** a0_list_i, sd_list** a1_list_i, int* jz_shell, int* l_shell);

void build_one_body_jumps_i_and_f(int n_s, int n_p, float mj_min, float mj_max, int num_mj, int n_sds_f, int n_sds_int, int* a1_array_f, wf_list** a0_list_i, sd_list** a1_list_i, int* jz_shell, int* l_shell);


#endif
