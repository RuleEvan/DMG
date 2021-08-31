#include <stdio.h>
#include "math.h"
#include "stdlib.h"
#include "stdint.h"
#include "string.h"
#include "time.h"
#include <gsl/gsl_sf.h>

#define HASH_SIZE 9781

// FILE SETUP
#define DENSITY_FILE "ne-mg_fermi_density"
#define WFN_FILE_INITIAL "ne20_basis.trwfn"
#define WFN_FILE_FINAL "mg20_basis.trwfn"
#define ORBIT_FILE "sd.sps"

#define MIN(a,b) ((a) < (b) ? (a):(b))
#define MAX(a,b) ((a) > (b) ? (a):(b))



