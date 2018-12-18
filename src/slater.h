#ifndef SLATER_H
#define SLATER_H
#include "angular.h"
unsigned int p_step(int n_s, int n_p, int *m_p);
unsigned int n_choose_k(int n, int k);
void orbitals_from_p(unsigned int p, int n_s, int n_p, int* orbitals);
void generate_single_particle_states();
unsigned int bin_from_p(int n_s, int n_p, unsigned int p);
float m_from_p(unsigned int p, int n_s, int n_p, int* m_shell);
unsigned int a_op(int n_s, int n_p, unsigned int p, int n_op, int* phase, int j_min);
unsigned int a_dag_op(int n_s, int n_p, unsigned int p, int n_op, int* phase);
unsigned int a_dag_a_op(int n_s, int n_p, unsigned int p, int n_a, int n_b, int* phase);
unsigned int a4_op(int n_s, int n_p, unsigned int p, int n_a, int n_b, int n_c, int n_d, int* phase);
unsigned int a4_op_b(int n_s, unsigned int b, int n_a, int n_b, int n_c, int n_d, int* phase);
unsigned int bin_phase_from_p(int n_s, int n_p, unsigned int p, int n_op, int* phase);
unsigned int bin_from_orbitals(int n_s, int n_p, int *orbitals);
void initialize_orbitals(int* n_shell, int* l_shell, int* j_shell, int* m_shell);
unsigned int count_set_bits(unsigned int n);
unsigned int p_from_binary(int n_s, int n_p, unsigned int b);
void orbitals_from_binary(int n_s, int n_p, unsigned int b, int *orbitals);
void generate_binomial_file();
unsigned int a2_op_b(int n_s, unsigned int b, int n_a, int n_b, int* phase);
int bit_count(unsigned int b);
int i_compare(const void * a, const void * b); 
unsigned int a_op_b(int n_s, unsigned int b, int n_a, int* phase);
int j_min_from_p(int n_s, int n_p, unsigned int p);
int max_mj(int n_s, int n_p, int *m_shell);
int min_mj(int n_s, int n_p, int *m_shell);
unsigned int get_num_sds(int n_s, int n_p);
#endif
