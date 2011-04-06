#include "Gromacs.h"
#include "SasAtom.h"
#include "SasAnalysis.h"

#include <string>
#include <iostream>
#include <cstring>

#include <boost/filesystem.hpp>

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
    gotTopology = false;
    gotTrajectory = false;
    cachedNFrames = 0;

    // Damn it! I can't handle errors raised inside this f*****g function,
    // because it simply crashes on a ERROR HANDLING FUNCTION, overriding
    // my SIGSEGV handler. Grrrr...    
    aps = gmx_atomprop_init();
  }

  Gromacs::~Gromacs()
  {
    operationThread.join();
#ifdef GMX45
    if(gotTrajectory)
      output_env_done(oenv);
#endif
    gmx_atomprop_destroy(aps);
  }
  
  void
  Gromacs::calculateSas()
  {
    operationThread = boost::thread(boost::ref(*this));
  }
  
  void
  Gromacs::operator ()()
  {
    bool bTop, bDGsol;
    real totarea, totvolume;
    int nsurfacedots;
    real *dgs_factor, *radius, *area, *surfacedots, dgsolv;
    atom_id* index;
    gmx_rmpbc_t gpbc;
    int nx;
    
    if(not (bTop = getTopology()))
      gmx_fatal(FARGS, "Could not read topology file.\n");
      
    //bDGsol = (strcmp(*(top.atoms.atomtype[0]),"?") != 0);
    // TODO: warnings about bDGsol
    // For now I really don't want anything related to DG solvatation!
    bDGsol = false;
    
    if(not getTrajectory())
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

    nx = grps->index[targetIndex + 1] - grps->index[targetIndex];
    index = new atom_id[nx];
    
    for(int i = 0; i < nx; i++)
      index[i] = grps->a[grps->index[targetIndex] + i];
    
    if(bDGsol)
      dgs_factor = new real[nx];
    radius = new real[natoms];
    
    for(int i = 0; i < natoms; i++)
    {
      gmx_atomprop_query( aps, epropVDW,
                          *(top.atoms.resinfo[top.atoms.atom[i].resind].name),
                          *(top.atoms.atomname[i]), radius + i);
      radius[i] += solSize;
    }
    
    if(bDGsol)
    {
      for(int i = 0; i < nx; i++)
      {
        int ii = index[i];
        if( !gmx_atomprop_query( aps, epropDGsol,
                            *(top.atoms.resinfo[top.atoms.atom[ii].resind].name),
                            *(top.atoms.atomname[ii]), dgs_factor + i))
          dgs_factor[i] = 0;
      }
    }
    
    gpbc = gmx_rmpbc_init(&top.idef, ePBC, natoms, box);
    SasAnalysis sasAnalysis(nx);
    
    do
    {
      gmx_rmpbc(gpbc, natoms, box, x);
      if(nsc_dclm_pbc(x, radius, nx, 24, FLAG_ATOM_AREA, &totarea,
                            &area, &totvolume, &surfacedots, &nsurfacedots,
                            index, ePBC, box) != 0)
        gmx_fatal(FARGS, "Something wrong in nsc_dclm_pbc");
        
      SasAtom atoms[nx];
      dgsolv = 0;
      for(int i = 0; i < nx; i++)
      {
        atoms[i].x = x[index[i]][0];
        atoms[i].y = x[index[i]][1];
        atoms[i].z = x[index[i]][2];
        atoms[i].sas = area[i];
        
        if(bDGsol)
          dgsolv += area[i]*dgs_factor[i];
      }
      
      sasAnalysis << atoms;
      
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

    if(bDGsol)
      delete[] dgs_factor;
    delete[] radius;
    delete[] gnames;
    delete[] grps->index;
    delete[] grps;
    
    //return true;
  }
  
  unsigned long
  Gromacs::getAtomsCount() const
  {
    if(gotTopology)
      return natoms;
    else
      return 0;
  }
  
  bool Gromacs::getTopology()
  {
    t_inputrec ir;
    matrix topbox;
    rvec *xtop;
    char title[1024];

#ifdef GMX45    
    read_tpx(tprName.c_str(), &ir, box, &natoms, 0, 0, 0, &mtop);
#else
    read_tpx( tprName.c_str(), &step, &t, &lambda, &ir, box, &natoms,
              0, 0, 0, &mtop);
#endif

    gotTopology = read_tps_conf( tprName.c_str(), title, &top, &ePBC, &xtop, 
                          NULL, topbox, FALSE);
    return gotTopology;
  }
  
  bool
  Gromacs::getTrajectory()
  {
#ifdef GMX45
    if(not gotTrajectory)
    {
      oenv = new output_env();
      output_env_init_default(oenv);
    }
    
    if((natoms = 
          read_first_x(oenv, &status, trjName.c_str(), &t, &x, box)) == 0)
#else
    if((natoms = 
          read_first_x(&status, trjName.c_str(), &t, &x, box)) == 0)
#endif
    {
      output_env_done(oenv);
      return gotTrajectory = false;
    }
    else
      return gotTrajectory = true;
  }
  
  unsigned int
  Gromacs::getFramesCount() const
  {
    if(cachedNFrames > 0)
      return cachedNFrames;
    
    namespace file = boost::filesystem;
    if(not file::exists(file::path(trjName)))
      return 0;
    
    unsigned int nFrames;
    output_env_t _oenv;
    t_trxstatus *_status;
    real _t;
    rvec* _x;
    matrix _box;
    int _natoms;
    
    _oenv = new output_env();
    output_env_init_default(_oenv);
    
    if((_natoms = read_first_x(_oenv,&_status, trjName.c_str(), &_t, &_x, 
                               _box)) == 0)
    {
      output_env_done(_oenv);
      return 0;
    }
    nFrames = 1;
    
    while(read_next_x(_oenv, _status, &_t, _natoms, _x, _box))
      nFrames++;
    
    output_env_done(_oenv);
    
    return cachedNFrames = nFrames;
  }
};
