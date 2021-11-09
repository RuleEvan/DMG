#include "density.h"

int main(int argc, char *argv[]) {
  clock_t start, end;
  double cpu_time;
  start = clock();
  if (argc != 2) {printf("Please supply only the parameter file name to the command line\n"); exit(0);}
  speedParams* sp = read_parameter_file(argv[1]);
  if ((sp->n_body == 1)) {
    one_body_density(sp);
  } else if ((sp->n_body == 2)) {
    two_body_density(sp);
  }
  end = clock();
  cpu_time = ((double) (end - start))/CLOCKS_PER_SEC;
  printf("Time: %g sec\n", cpu_time);
  return 0;
}
