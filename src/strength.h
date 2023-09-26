#ifndef STRENGTH_H
#define STRENGTH_H
#include "density.h"

 
void two_body_pivot(speedParams* sp);



/*   -------------------------
 *   Two-body routines       
 *   -------------------------
 */

// Two-body node trace routines

void trace_two_body_nodes_dmtp2_strength(int a, int b, int c, int d, int num_mj_p_i, float mj_min_p_i, int num_mj_n_i, float mj_min_n_i, sd_list** a2_list_i, sd_list** a2_list_f, wfnData* wd, int i_op, eigen_list* transition, double* density);

double compute_matrix_element_F0(int in1p, int ij1p, int in2p, int ij2p, int ij12p, int in1, int ij1, int in2, int ij2, int ij12, int it12); 
int get_l(int n, int j);

double double_fermi(int a, int b, int c, int d, int* n_list, int* l_list, int* j_list, int* mj_list); 
double double_gamow_teller(int a, int b, int c, int d, int* n_list, int* l_list, int* j_list, int* mj_list); 
double double_gamow_teller2(int a, int b, int c, int d, int* n_list, int* l_list, int* j_list, int* mj_list); 

double three_j(double j1, double j2, double j3, double m1, double m2, double m3);
double six_j(double j1, double j2, double j3, double j4, double j5, double j6);
double nine_j(double j11, double j12, double j13, double j21, double j22, double j23, double j31, double j32, double j33);
double triangle(double a, double b, double c);

// Two-body jump routines

void build_two_body_jumps_dmt0(int n_s, int n_p, float mj_min, float mj_max, int num_mj, int n_sds_i, int n_sds_int1, int n_sds_int2, int*a1_array_f, int* a2_array_f, wf_list** a0_list_i, sd_list** a1_list_i, sd_list** a2_list_i, int* jz_shell, int* l_shell);

void build_two_body_jumps_dmtp1(int n_s, int n_proton_i, int num_mj_p_i, float mj_min_p_i, float mj_max_p_i, int n_proton_f, int num_mj_p_f, float mj_min_p_f, float mj_max_p_f, sd_list** p1_list_i, sd_list** p1_list_f, int* p2_array_f, int n_neutron_i, int num_mj_n_i, float mj_min_n_i, float mj_max_n_i, int n_neutron_f, int num_mj_n_f, float mj_min_n_f, float mj_max_n_f, sd_list** n1_list_i, sd_list** n2_list_i, int* n1_array_f, int* jz_shell, int* l_shell);

void build_two_body_jumps_dmtp2(int n_s, int n_proton, int num_mj_p_i, float mj_min_p_i, float mj_max_p_i, int num_mj_p_f, float mj_min_p_f, float mj_max_p_f, int n_sds_p_f, sd_list** p2_list_f, int n_neutron, int num_mj_n_i, float mj_min_n_i, float mj_max_n_i, int num_mj_n_f, float mj_min_n_f, float mj_max_n_f, int n_sds_n_i, sd_list** n2_list_i, int* jz_shell, int* l_shell);

double compute_matrix_element_GT0(int in1p, int ij1p, int in2p, int ij2p, int ij12p, int in1, int ij1, int in2, int ij2, int ij12, int it12);
double compute_matrix_element_GT2(int in1p, int ij1p, int in2p, int ij2p, int ij12p, int in1, int ij1, int in2, int ij2, int ij12, int it12);

double gamow_teller(int a, int b, int* n_list, int* l_list, int* j_list, int* mj_list); 
void trace_one_body_nodes_dmtp1_strength(int a, int b, int num_mj_p_i, float mj_min_p_i, int num_mj_n_i, float mj_min_n_i, sd_list** n1_list_i, sd_list** p1_list_f, wfnData* wd, int i_op, eigen_list *transition, double* density); 
void one_body_pivot(speedParams* sp); 
#endif
