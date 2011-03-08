#include "Gromacs.h"
#include <string>
#include <iostream>

using namespace std;

namespace Gromacs
{
  Gromacs::Gromacs()
  {
    init();
  }

  Gromacs::Gromacs( const string& trajectoryFileName,
                    const string& topologyFileName)
  {
    trjName = trajectoryFileName;
    tprName = topologyFileName;
    init();
  }

  void
  Gromacs::init()
  {
#ifndef GMX45
    // FIXME: Test versions until 4.5
    static const char* argv[] = { "gromacs" };
    static int argc = 1;
    cr = init_par(&argc, &argv);
#endif
    
    aps = gmx_atomprop_init();
  }

  Gromacs::~Gromacs()
  {
    gmx_atomprop_destroy(aps);
  }
  
  bool
  Gromacs::calculateSas() const
  {

#ifdef GMX45
    output_env_t oenv;
    t_trxstatus *status;
#else
    t_commrec *cr;
    int status, step;
    real t, lambda;
#endif

    int natoms;
    t_inputrec ir;
    matrix box;
    gmx_mtop_t mtop;
    t_trxframe fr;
    
#ifdef GMX45
    oenv = new output_env;
    output_env_init_default(oenv);
    read_tpx(tprName.c_str(), &ir, box, &natoms, 0, 0, 0, &mtop);
#else
    read_tpx( tprName.c_str(), &step, &t, &lambda, &ir, box, &natoms,
              0, 0, 0, &mtop);
#endif
    
#ifdef GMX45
    if(!read_first_frame(oenv, &status, trjName.c_str(), &fr, TRX_NEED_X))
#else
    if(!read_first_frame(&status, trjName.c_str(), &fr, TRX_NEED_X))
#endif
      gmx_fatal(FARGS, "The trajectory file is invalid.\n");
    
    do
    {
      cout <<  endl << fr.natoms << " atoms read." << endl;
#ifdef GMX45
    } while(read_next_frame(oenv, status, &fr));
#else
    } while(read_next_frame(status, &fr));
#endif

    close_trx(status);
#ifdef GMX45    
    output_env_done(oenv);
#endif

    return true;
  }
};
