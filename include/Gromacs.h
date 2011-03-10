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
#include <rmpbc.h>
};

/* Taken from src/tools/nsc.h, Gromacs 4.5.2 */
#define FLAG_ATOM_AREA  04
extern "C" int nsc_dclm_pbc(rvec *coords, real *radius, int nat,
                        int  densit, int mode,
                        real *value_of_area, real **at_area,
                        real *value_of_vol,
                        real **lidots, int *nu_dots,
                        atom_id index[],int ePBC,matrix box);

using namespace std;

namespace Gromacs
{
  class Gromacs
  {
  public:
    Gromacs(float solventSize = 0.14);
    ~Gromacs();
    Gromacs(const string& trajectoryFileName,
            const string& topologyFileName,
            float solventSize = 0.14);
    // FIXME: It will have to return an object Molecular Dynamics with SAS
    // FIXME: additional informations
    bool calculateSas() const;

  private:
    gmx_atomprop_t aps;
    string trjName, tprName;
    float solSize;
    
    void init(float solventSize);
  };
};

#endif

