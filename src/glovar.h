#include <stdio.h>
#include "math.h"
#include "stdlib.h"
#include "stdint.h"
#include "string.h"
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sf_expint.h>
#include <gsl/gsl_sf_result.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_sf_bessel.h>

// Isotope setup
#define A_NUC 76 // Atomic Mass
#define NUM_SHELLS 259
#define A_FACTOR 9.155 // [MeV] Average nuclear excitation energy
#define B_OSC 0.8
#define Z_ATOM 32 // Atomic Number
#define HASH_SIZE 9781

// FILE SETUP
#define DENSITY_FILE "ne-mg_fermi_density"
#define WFN_FILE_INITIAL "ne20_basis.trwfn"
#define WFN_FILE_FINAL "mg20_basis.trwfn"
#define ORBIT_FILE "sd.sps"

// Technical parameters
#define R_SPL_MIN 0.0001
#define R_SPL_MAX 100.0
#define R_SPL_STEPS 1000
#define COR_FAC 1
#define POTENTIAL 1

// Physical constants
#define HIGGS_VEV 246.0 // [GeV]
#define V_UD 0.97427
#define ALPHA_FS 0.007297352
#define R_NUC 1.2 // [fm]
#define M_NUC 2.19 // nuclear matrix element
#define M_ELECTRON 0.511 // [MeV]
#define M_NEUTRON 939.57 // [MeV]
#define E_BETA 2.55 // [MeV]
#define G_AXIAL 1.275 // Axial coupling
#define G_MAGNETIC 4.7
#define G_TENSOR 0.99
#define G_VECTOR 1.0
#define E_HADRON 1 // Hadronization scale [GeV] 
#define DELTA_SUPP 1.0/30.0 // Suppression factor
#define MU 0.6995 // [fm]^-1
#define KAPPA_1 3.7
#define LAMBDA_A 1040.0 //[MeV]
#define LAMBDA_V 850.0 // [MeV]
#define PION_MASS 139.57 // [MeV] Charged pion mass

#define MIN(a,b) ((a) < (b) ? (a):(b))
#define MAX(a,b) ((a) > (b) ? (a):(b))



