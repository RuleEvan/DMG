#include <stdio.h>
#include "math.h"
#include "stdlib.h"
#include "stdint.h"
#include "string.h"
#include "time.h"
#include <gsl/gsl_sf.h>

// Set the size of the hash used for storing wavefunction coefficients
#define HASH_SIZE 9781

// Choose print format for one-body density matrices
// 0 is GLASGOW format
// 1 is BIGSTICK format
#define I_FORMAT 0

// FILE SETUP

#define MIN(a,b) ((a) < (b) ? (a):(b))
#define MAX(a,b) ((a) > (b) ? (a):(b))



