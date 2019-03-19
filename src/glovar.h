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

#define HASH_SIZE 9781

// FILE SETUP
#define DENSITY_FILE "ne-mg_fermi_density"
#define WFN_FILE_INITIAL "ne20_basis.trwfn"
#define WFN_FILE_FINAL "mg20_basis.trwfn"
#define ORBIT_FILE "sd.sps"

#define MIN(a,b) ((a) < (b) ? (a):(b))
#define MAX(a,b) ((a) > (b) ? (a):(b))



