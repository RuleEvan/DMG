#ifndef SLATER_H
#define SLATER_H
#include "angular.h"
unsigned int p_step(int n_s, int n_p, int *m_p);
unsigned long n_choose_k(unsigned int n, unsigned int k);
//void orbitals_from_p(unsigned int p, int n_s, int n_p, int* orbitals);
float m_from_p(unsigned int p, int n_s, int n_p, int* m_shell);
int parity_from_p(unsigned int p, int n_s, int n_p, int* l_shell);
unsigned int a_op(int n_s, int n_p, unsigned int p, int n_op, int* phase, int j_min);
unsigned int a_op_dag(int n_s, int n_p, unsigned int p, int n_op, int* phase, int j_min);
void generate_binomial_file();
int j_min_from_p(int n_s, int n_p, unsigned int p);
float max_mj(int n_s, int n_p, int *m_shell);
float min_mj(int n_s, int n_p, int *m_shell);
unsigned int get_num_sds(int n_s, int n_p);
int w_from_p(unsigned int p, int n_s, int n_p, int* w_shell);
float min_w(int n_s, int n_p, int *w_shell);

int test_n_choose_k();
int test_slater_routines();

#endif
