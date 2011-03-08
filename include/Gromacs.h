#ifndef _GROMACS_H
#define _GROMACS_H

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <string>

extern "C" {
#include <atomprop.h>
#include <statutil.h>
#include <tpxio.h>
#include <mtop_util.h>
#include <main.h>
};

using namespace std;

namespace Gromacs
{
  class Gromacs
  {
  public:
    Gromacs();
    ~Gromacs();
    Gromacs(const string& trajectoryFileName, const string& topologyFileName);
    // FIXME: It will have to return an object Molecular Dynamics with SAS
    // FIXME: additional informations
    bool calculateSas() const;

  private:
    gmx_atomprop_t aps;
    string trjName, tprName;
    
    void init();
  };
};

#endif

