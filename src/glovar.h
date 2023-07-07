#include "stdio.h"
#include "math.h"
#include "stdlib.h"
#include "stdint.h"
#include "string.h"
#include "time.h"
#include <gsl/gsl_sf.h>


// Set the size of the hash used for storing wavefunction coefficients
#define HASH_SIZE 62273
// Choose print format for one-body density matrices
// 0 is GLASGOW format
// 1 is BIGSTICK format
#define I_FORMAT 0
#define I_VERBOSE 1

// Macros for MIN and MAX functions
#define MIN(a,b) ((a) < (b) ? (a):(b))
#define MAX(a,b) ((a) > (b) ? (a):(b))



