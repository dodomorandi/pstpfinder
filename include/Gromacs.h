#ifndef _GROMACS_H
#define _GROMACS_H

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

extern "C" {
#include <atomprop.h>
#include <statutil.h>
}

namespace Gromacs
{
  class Gromacs
  {
  public:
    Gromacs();
    ~Gromacs();
  private:
    #ifdef GMX45
    output_env_t oenv;
    t_trxstatus *status;
    #else
    t_commrec *cr;
    int status, step;
    real t, lambda;
    #endif
    
    gmx_atomprop_t aps;
  };
};

#endif

