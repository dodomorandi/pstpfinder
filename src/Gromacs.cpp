#include "Gromacs.h"
#include <string>
#include <iostream>
#include <cstring>

using namespace std;

namespace Gromacs
{
  Gromacs::Gromacs(float solventSize)
  {
    init(solventSize);
  }

  Gromacs::Gromacs( const string& trajectoryFileName,
                    const string& topologyFileName, 
                    float solventSize)
  {
    trjName = trajectoryFileName;
    tprName = topologyFileName;
    init(solventSize);
  }

  void
  Gromacs::init(float solventSize)
  {
#ifndef GMX45
    // FIXME: Test versions until 4.5
    static const char* argv[] = { "gromacs" };
    static int argc = 1;
    cr = init_par(&argc, &argv);
#endif
    
    solSize = solventSize;
    sasTarget = "Protein";
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

    bool bTop, bDGsol;
    bool *bOut;
    t_inputrec ir;
    matrix box;
    gmx_mtop_t mtop;
    matrix topbox;
    rvec *xtop, *x;
    real t, totarea, totvolume, tarea;
    int natoms, ePBC, nsurfacedots;
    t_topology top;
    char title[1024];
    real *dgs_factor, *radius, *area, *surfacedots, dgsolv;
    atom_id* index[2];
    gmx_rmpbc_t gpbc;
    int nx[2];
    
#ifdef GMX45
    oenv = new output_env();
    output_env_init_default(oenv);
    read_tpx(tprName.c_str(), &ir, box, &natoms, 0, 0, 0, &mtop);
#else
    read_tpx( tprName.c_str(), &step, &t, &lambda, &ir, box, &natoms,
              0, 0, 0, &mtop);
#endif
    
    bTop = read_tps_conf( tprName.c_str(), title, &top, &ePBC, &xtop, 
                          NULL, topbox, FALSE);
    
    bDGsol = (strcmp(*(top.atoms.atomtype[0]),"?") != 0);
    // TODO: warnings about bDGsol
    
#ifdef GMX45
    if((natoms = 
          read_first_x(oenv, &status, trjName.c_str(), &t, &x, box)) == 0)
#else
    if((natoms = 
          read_first_x(&status, trjName.c_str(), &t, &x, box)) == 0)
#endif
      gmx_fatal(FARGS, "Could not read coordinates from statusfile.\n");
    
    // TODO: index file handling and integration
    
    // The original code was simply:
    // get_index(&top.atoms, 0, 2, nx, index, grpname);
    // but now I need an automatic selection or a GUI.
    // For this purpose it's necessary to partially use index.c code in GROMACS
    
    // These will be reallocated by Gromacs libraries.
    // They must be allocated on the heap.
    // I think it'g good to allocate them as arrays, because they will probably
    // be reallocated. It would be easier to debug later.
    t_blocka* grps;
    grps = new t_blocka[1]();
    grps->index = new atom_id[1]();
    char*** gnames;
    gnames = new char**[1]();
    int targetIndex = -1;
    
    analyse(&top.atoms, grps, gnames, FALSE, FALSE);
    
    for(char** i = *gnames; i < *gnames + grps->nr; i++)
    {
      if(*i == sasTarget)
      {
        targetIndex = (int)(i - *gnames);
        break;
      }
    }
    // TODO: Return in case of missing target index

    nx[0] = nx[1] = grps->index[targetIndex + 1] - grps->index[targetIndex];
    index[0] = new atom_id[nx[0]];
    index[1] = new atom_id[nx[1]];
    
    for(int i = 0; i < nx[0]; i++)
      index[0][i] = index[1][i] = grps->a[grps->index[targetIndex] + i];
    
    bOut = new bool[natoms];
    for(int* i = index[1]; i < index[1] + nx[1]; i++)
      bOut[*i] = TRUE;
    
    dgs_factor = new real[nx[0]];
    radius = new real[natoms];
    
    for(int i = 0; i < natoms; i++)
    {
      gmx_atomprop_query( aps, epropVDW,
                          *(top.atoms.resinfo[top.atoms.atom[i].resind].name),
                          *(top.atoms.atomname[i]), radius + i);
      radius[i] += solSize;
    }
    for(int i = 0; i < nx[0]; i++)
    {
      int ii = index[0][i];
      if( !gmx_atomprop_query( aps, epropVDW,
                          *(top.atoms.resinfo[top.atoms.atom[ii].resind].name),
                          *(top.atoms.atomname[ii]), dgs_factor + i))
        dgs_factor[i] = 0;
    }
    
    gpbc = gmx_rmpbc_init(&top.idef, ePBC, natoms, box);
    
    do
    {
      gmx_rmpbc(gpbc, natoms, box, x);
      if(nsc_dclm_pbc(x, radius, nx[0], 24, FLAG_ATOM_AREA, &totarea,
                            &area, &totvolume, &surfacedots, &nsurfacedots,
                            index[0], ePBC, box) != 0)
        gmx_fatal(FARGS, "Something wrong in nsc_dclm_pbc");
        
        tarea = 0;
        dgsolv = 0;
        for(int i = 0; i < nx[0]; i++)
        {
          if(bOut[index[0][i]])
          {
            tarea += area[i];
            if(bDGsol)
              dgsolv += area[i]*dgs_factor[i];
          }
        }
      cout << "At " << t << " ps area is " << tarea << endl;
      
      delete[] area;
      if(nsurfacedots > 0)
        delete[] surfacedots;
    }  
#ifdef GMX45
    while(read_next_x(oenv, status, &t, natoms, x, box));
#else
    while(read_next_x(status, &t, natoms, x, box));
#endif

    gmx_rmpbc_done(gpbc);

    close_trx(status);
#ifdef GMX45    
    output_env_done(oenv);
#endif

    delete[] bOut;
    delete[] dgs_factor;
    delete[] radius;
    delete[] gnames;
    delete[] grps->index;
    delete[] grps;
    
    return true;
  }
};
