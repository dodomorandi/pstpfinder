#include "Gromacs.h"

Gromacs::Gromacs::Gromacs()
{
  #ifdef GMX45
  oenv = new output_env;
  #else
  // FIXME: Test versions until 4.5
  static char* argv[] = { "gromacs" };
  static int argc = 1;
  cr = init_par(&argc, &argv);
  #endif
  
  aps = gmx_atomprop_init();
}

Gromacs::Gromacs::~Gromacs()
{
  gmx_atomprop_destroy(aps);
  #ifdef GMX45
  output_env_done(oenv);
  #endif
}
