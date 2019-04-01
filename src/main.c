#include "density.h"

int main(int argc, char *argv[]) {
  if (argc != 2) {printf("Please supply only the parameter file name to the command line\n"); exit(0);}
  speedParams* sp = read_parameter_file(argv[1]);
  if ((sp->n_body == 1) && (sp->spec_dep == 0)) {
    one_body_density(sp);
  } else if ((sp->n_body == 1) && (sp->spec_dep == 1)) {
    one_body_density_spec(sp);
  } else if ((sp->n_body == 2) && (sp->spec_dep == 0)) {
    two_body_density(sp);
  } else if ((sp->n_body == 2) && (sp->spec_dep == 1)) {
    two_body_density_spec(sp);
  }
  return 0;
}
