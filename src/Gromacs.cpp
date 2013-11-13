/*
 *  This file is part of PSTP-finder, an user friendly tool to analyze GROMACS
 *  molecular dynamics and find transient pockets on the surface of proteins.
 *  Copyright (C) 2011 Edoardo Morandi.
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "Gromacs.h"
#include "SasAtom.h"
#include "SasAnalysis.h"
#include "utils.h"
#include "Protein.h"

#include <string>
#include <iostream>
#include <fstream>
#include <cstring>
#include <mutex>

#include <gromacs/tpxio.h>
#include <gromacs/mtop_util.h>
#include <gromacs/main.h>
#include <gromacs/rmpbc.h>
#include <gromacs/xtcio.h>
#include <gromacs/typedefs.h>
#include <gromacs/princ.h>
#include <gromacs/do_fit.h>
#include <gromacs/smalloc.h>
#include <gromacs/vec.h>
#include <gromacs/statutil.h>
#include <gromacs/xdrf.h>
#include <gromacs/confio.h>

namespace PstpFinder
{
  static std::mutex nsc_dclm_pbc_mutex;

  Gromacs::Gromacs(float solventSize)
  {
    init(solventSize);
  }

  Gromacs::Gromacs(const std::string& trajectoryFileName,
                   const std::string& topologyFileName, float solventSize)
  {
    trjName = trajectoryFileName;
    tprName = topologyFileName;
    init(solventSize);
  }

  Gromacs::Gromacs(const Gromacs& gromacs)
  {
    init(gromacs.solSize);

    trjName = gromacs.trjName;
    tprName = gromacs.tprName;
    sasTarget = gromacs.sasTarget;
    _usePBC = gromacs._usePBC;

    cachedNFrames = gromacs.cachedNFrames;
    averageStructure = gromacs.averageStructure;
    _begin = gromacs._begin;
    _end = gromacs._end;
    timeStepCached = gromacs.timeStepCached;
  }

  void
  Gromacs::init(float solventSize)
  {
#ifndef GMX45
    // FIXME: Test versions until 4.5
    static const char* argv[] =
    { "gromacs"};
    static int argc = 1;
    cr = init_par(&argc, &argv);
#endif

    solSize = solventSize;
    sasTarget = "Protein";
    gotTopology = false;
    gotTrajectory = false;
    readyToGetX = false;
    currentFrame = 0;
    cachedNFrames = 0;
    _begin = -1;
    _end = -1;
    timeStepCached = 0;
    abortFlag = false;
    _usePBC = true;

    // Damn it! I can't handle errors raised inside this f*****g function,
    // because it simply crashes on a ERROR HANDLING FUNCTION, overriding
    // my SIGSEGV handler. Grrrr...    
    aps = gmx_atomprop_init();
  }

  Gromacs::~Gromacs()
  {
    if(operationThread.joinable())
      operationThread.join();
#ifdef GMX45
    if(gotTrajectory)
    {
      output_env_done(oenv);
      close_trx(status);
    }

    if(gotTopology)
      sfree(xtop);
#endif
    gmx_atomprop_destroy(aps);
  }

  template<typename Stream>
  void
  Gromacs::calculateSas(Session<Stream>& session)
  {
    operationThread = std::thread(
        bind(&Gromacs::__calculateSas<Stream>, std::ref(*this), std::ref(session)));
  }

  void
  Gromacs::calculateAverageStructure()
  {
    operationThread = std::thread(std::bind(&Gromacs::__calculateAverageStructure,
                                  std::ref(*this)));
  }

  void
  Gromacs::waitOperation()
  {
    if(operationThread.joinable())
      operationThread.join();
  }

  template<typename Stream>
  void
  Gromacs::__calculateSas(Session<Stream>& session)
  {
    bool bDGsol;
    real totarea, totvolume;
    int nsurfacedots;
    real *dgs_factor, *radius, *area = 0, *surfacedots = 0, dgsolv;
    std::vector<atom_id> index;
    gmx_rmpbc_t gpbc = nullptr;
    int nx;

    if(not gotTopology and not getTopology())
      gmx_fatal(FARGS, "Could not read topology file.\n");

    /* Just a precaution -- we could change this trough SasAnalysis */
    operationMutex.lock();
    currentFrame = 0;
    operationMutex.unlock();

    /* bDGsol = (strcmp(*(top.atoms.atomtype[0]),"?") != 0); */
    /* TODO: warnings about bDGsol */
    /* For now I really don't want anything related to DG solvatation! */
    bDGsol = false;

    if(not gotTrajectory and not getTrajectory())
      gmx_fatal(FARGS, "Could not read coordinates from statusfile.\n");

    /* TODO: index file handling and integration */

    index = getGroup("Protein");
    nx = index.size();

    if(bDGsol)
      dgs_factor = new real[nx];
    radius = new real[natoms];

    for(int i = 0; i < natoms; i++)
    {
      gmx_atomprop_query(aps, epropVDW,
                         *(top.atoms.resinfo[top.atoms.atom[i].resind].name),
                         *(top.atoms.atomname[i]), radius + i);
      radius[i] += solSize;
    }

    if(bDGsol)
    {
      for(int i = 0; i < nx; i++)
      {
        int ii = index[i];
        if(!gmx_atomprop_query(
            aps, epropDGsol,
            *(top.atoms.resinfo[top.atoms.atom[ii].resind].name),
            *(top.atoms.atomname[ii]), dgs_factor + i))
          dgs_factor[i] = 0;
      }
    }

    if(_usePBC)
      gpbc = gmx_rmpbc_init(&top.idef, ePBC, natoms, fr.box);

    SasAnalysis<Stream> sasAnalysis(nx, *this, session);
    {
      unsigned int readFrames(sasAnalysis.getReadFrames());
      while(readFrames > currentFrame)
      {
        //gmx_rmpbc(gpbc, natoms, box, x);
        operationMutex.lock();
        currentFrame++;
        wakeCondition.notify_all();
        operationMutex.unlock();
        if(not readNextX())
          break;
      }
    }

    do
    {
      if(abortFlag)
      {
        session.abort();
        break;
      }
      if(_usePBC)
        gmx_rmpbc(gpbc, natoms, fr.box, fr.x);

      int nsc_dclm_pdc_result;
      {
        std::unique_lock<std::mutex> lock(nsc_dclm_pbc_mutex);
        nsc_dclm_pdc_result = nsc_dclm_pbc(fr.x, radius, nx, 24, FLAG_ATOM_AREA,
                                           &totarea, &area, &totvolume,
                                           &surfacedots, &nsurfacedots,
                                           index.data(), ePBC,
                                           _usePBC ? fr.box : nullptr);
      }
      if(nsc_dclm_pdc_result != 0)
        gmx_fatal(FARGS, "Something wrong in nsc_dclm_pbc");

      SasAtom atoms[nx];
      dgsolv = 0;
      for(int i = 0; i < nx; i++)
      {
        atoms[i].x = fr.x[index[i]][0];
        atoms[i].y = fr.x[index[i]][1];
        atoms[i].z = fr.x[index[i]][2];
        atoms[i].sas = area[i];

        if(bDGsol)
          dgsolv += area[i] * dgs_factor[i];
      }
      if(abortFlag)
      {
        session.abort();
        break;
      }

      operationMutex.lock();
      sasAnalysis << atoms;
      currentFrame++;
      wakeCondition.notify_all();
      operationMutex.unlock();

      if(area)
      {
        sfree(area);
        area = NULL;
      }
      if(surfacedots)
      {
        sfree(surfacedots);
        surfacedots = NULL;
      }
    }
    while(readNextX());

    if(_usePBC)
      gmx_rmpbc_done(gpbc);

    output_env_done(oenv);
    close_trx(status);
    gotTrajectory = false;

    if(bDGsol)
      delete[] dgs_factor;
    delete[] radius;
  }

  const Protein<>&
  Gromacs::__calculateAverageStructure()
  {
    std::vector<atom_id> index;
    int npdbatoms, isize, count;
    rvec xcm;
    matrix pdbbox;
    real *w_rls;
    real invcount;
    double *xav, *rmsf;
    double** U;
    gmx_rmpbc_t gpbc = nullptr;
    int statusCount = 0;
    averageStructure = Protein<>();

    if(not getTopology())
      gmx_fatal(FARGS, "Could not read topology file.\n");

    if(not gotTrajectory and not getTrajectory())
      gmx_fatal(FARGS, "Could not read coordinates from statusfile.\n");

    currentFrame = 0;

    index = getGroup("Protein");
    isize = index.size();

    w_rls = new real[top.atoms.nr];
    for(std::vector<atom_id>::const_iterator i = index.begin(); i < index.end(); i++)
      w_rls[*i] = top.atoms.atom[*i].m;

    xav = new double[isize * DIM]();
    U = new double*[isize];
    for(double** i = U; i < U + isize; i++)
      *i = new double[DIM * DIM]();
    rmsf = new double[isize];

    npdbatoms = top.atoms.nr;
    snew(top.atoms.pdbinfo, npdbatoms);
    copy_mat(fr.box, pdbbox);

    sub_xcm(xtop, isize, index.data(), top.atoms.atom, xcm, FALSE);
    if(_usePBC)
      gpbc = gmx_rmpbc_init(&top.idef, ePBC, natoms, fr.box);

    count = 0;
    do
    {
      if(abortFlag)
        break;
      if(_usePBC)
        gmx_rmpbc(gpbc, natoms, fr.box, fr.x);
      sub_xcm(fr.x, isize, index.data(), top.atoms.atom, xcm, FALSE);
      do_fit(natoms, w_rls, xtop, fr.x);

      if(abortFlag)
        break;
      for(int i = 0; i < isize; i++)
      {
        atom_id aid = index[i];
        for(int d = 0; d < DIM; d++)
        {
          xav[i * DIM + d] += fr.x[aid][d];
          for(int m = 0; m < DIM; m++)
            U[i][d * DIM + m] += fr.x[aid][d] * fr.x[aid][m];
        }
      }

      count++;
      operationMutex.lock();
      currentFrame = ++statusCount;
      wakeCondition.notify_all();
      operationMutex.unlock();
    }
    while(readNextX());

    if(abortFlag)
      return averageStructure;

    if(_usePBC)
      gmx_rmpbc_done(gpbc);
    invcount = 1.0 / count;
    for(int i = 0; i < isize; i++)
    {
      for(int d = 0; d < DIM; d++)
        xav[i * DIM + d] *= invcount;

      for(int d = 0; d < DIM; d++)
        for(int m = 0; m < DIM; m++)
          U[i][d * DIM + m] = U[i][d * DIM + m] * invcount
                              - xav[i * DIM + d] * xav[i * DIM + m];
    }

    for(int i = 0; i < isize; i++)
      rmsf[i] = U[i][XX * DIM + XX] + U[i][YY * DIM + YY] + U[i][ZZ * DIM + ZZ];

    for(double** i = U; i < U + isize; i++)
      delete[] *i;
    delete[] U;

    if(abortFlag)
      return averageStructure;
    for(int i = 0; i < isize; i++)
      top.atoms.pdbinfo[index[i]].bfac = 800 * M_PI * M_PI / 3.0 * rmsf[i];

    Residue<> res;
    for(int i = 0; i < isize; i++)
    {
      if(abortFlag)
        return averageStructure;
      ProteinAtom atom;
      atom.index = i + 1;
      atom.setAtomType(*top.atoms.atomname[index[i]], false);

      atom.x = xcm[0] + xav[i * DIM];
      atom.y = xcm[1] + xav[i * DIM + 1];
      atom.z = xcm[2] + xav[i * DIM + 2];

      /*
      if(top.atoms.pdbinfo)
      {
        atom.bFactor = top.atoms.pdbinfo[index[i]].bfac;
        atom.occupancy = top.atoms.pdbinfo[index[i]].occup;
      }
      else
      {
        atom.bFactor = 0.0;
        atom.occupancy = 1.0;
      }
      */

      int resind = top.atoms.atom[index[i]].resind;

      if(res.atoms.size() == 0 or top.atoms.resinfo[resind].nr != res.index)
      {
        if(res.atoms.size() != 0)
        {
          averageStructure.appendResidue(res);
          res = Residue<>();
        }

        res.index = top.atoms.resinfo[resind].nr;
        res.type = Residue<>::getTypeByName(
            *top.atoms.resinfo[resind].name);
        res.chain = 'A';
      }
      res.atoms.push_back(std::move(atom));

      operationMutex.lock();
      currentFrame = (float) getFramesCount() / isize * i;
      wakeCondition.notify_all();
      operationMutex.unlock();
    }
    averageStructure.appendResidue(res);

    delete[] rmsf;
    delete[] xav;
    delete[] w_rls;

    return averageStructure;
  }

  const Protein<>&
  Gromacs::getAverageStructure() const
  {
    return averageStructure;
  }

  void
  Gromacs::setAverageStructure(Protein<> structure)
  {
    averageStructure = structure;
  }

  unsigned long
  Gromacs::getAtomsCount() const
  {
    if(gotTopology)
      return natoms;
    else
      return 0;
  }

  std::vector<atom_id>
  Gromacs::getGroup(const std::string& groupName)
  {
    if(not gotTopology and not getTopology())
      abort();

    return const_cast<const Gromacs&>(*this).getGroup(groupName);
  }

  std::vector<atom_id>
  Gromacs::getGroup(const std::string& groupName) const
  {
    // The original code was simply:
    // get_index(&top.atoms, 0, 2, nx, index, grpname);
    // but now I need an automatic selection or a GUI.
    // For this purpose it's necessary to partially use index.c code in GROMACS

    if(not gotTopology)
      throw;

    std::vector<atom_id> group;
    int size;

    t_blocka* grps;
    snew(grps, 1);
    snew(grps->index, 1);
    char*** gnames;
    snew(gnames, 1);
    int targetIndex = -1;

    t_atoms* m_atoms = new t_atoms;
    memcpy(m_atoms, &top.atoms, sizeof(t_atoms));

    analyse(m_atoms, grps, gnames, FALSE, FALSE);

    for(char** i = *gnames; i < *gnames + grps->nr; i++)
    {
      if(groupName.compare(*i) == 0)
      {
        targetIndex = (int) (i - *gnames);
        break;
      }
    }
    // TODO: Return in case of missing target index

    size = grps->index[targetIndex + 1] - grps->index[targetIndex];
    group.reserve(size);

    for(int i = 0; i < size; i++)
      group.push_back(grps->a[grps->index[targetIndex] + i]);

    sfree(gnames);
    sfree(grps->index);
    sfree(grps->a);
    sfree(grps);
    delete m_atoms;

    return group;
  }

  bool
  Gromacs::getTopology()
  {
    t_inputrec ir;
    matrix topbox;
    char title[1024];
    std::string extension = file_extension(tprName);
    std::transform(std::begin(extension), std::end(extension), std::begin(extension), ::tolower);

    if(extension == ".tpr")
    {
#ifdef GMX45
    read_tpx(tprName.c_str(), &ir, fr.box, &natoms, 0, 0, 0, &mtop);
#else
    read_tpx( tprName.c_str(), &step, &t, &lambda, &ir, box, &natoms,
        0, 0, 0, &mtop);
#endif
    }
    else if(extension == ".pdb")
    {
      get_stx_coordnum(tprName.c_str(), &natoms);
    }
    else
    {
      gotTopology = false;
      return gotTopology;
    }

    read_tps_conf(tprName.c_str(), title, &top, &ePBC, &xtop,
                                NULL, topbox, FALSE);

    gotTopology = true;
    return gotTopology;
  }

  bool
  Gromacs::getTrajectory()
  {
#ifdef GMX45
    if(not gotTrajectory)
    {
      snew(oenv, 1);
      output_env_init_default(oenv);
    }
    else
      close_trx(status);

    cachedNFrames = 0;
    timeStepCached = 0;
    currentFrame = 0;

    if(not read_first_frame(oenv, &status, trjName.c_str(), &fr,
                                  TRX_NEED_X))
#else
    if((natoms =
            read_first_x(&status, trjName.c_str(), &t, &x, box)) == 0)
#endif
    {
      output_env_done(oenv);
      readyToGetX = false;
      return gotTrajectory = false;
    }
    else
    {
#ifdef GMX45
      natoms = fr.natoms;
#endif
      readyToGetX = true;
      return gotTrajectory = true;
    }
  }

  std::string
  Gromacs::getTrajectoryFile() const
  {
    return trjName;
  }

  std::string
  Gromacs::getTopologyFile() const
  {
    return tprName;
  }

  bool
  Gromacs::readNextX()
  {
    bool out;

    if(not gotTrajectory or not readyToGetX)
      if(not getTrajectory())
        throw;

#ifdef GMX45
    out = read_next_frame(oenv, status, &fr);
#else
    out = read_next_x(status, &t, natoms, x, box);
#endif

    if(not out)
      readyToGetX = false;

    return out;
  }

  unsigned int
  Gromacs::getFramesCount() const
  {
    if(cachedNFrames > 0)
      return cachedNFrames;
    else if(_begin != -1 and _end != -1)
    {
      cachedNFrames = (_end - _begin) / getTimeStep();
      return cachedNFrames;
    }

    if(not exists(trjName))
      return 0;

    int nFrames = 0;
    std::string extension = file_extension(trjName);
    transform(std::begin(extension), std::end(extension), std::begin(extension), ::tolower);
    if(extension == "pdb")
    {
      /* I use old std legacy method for pdb files */
      output_env_t _oenv;
      t_trxstatus *_status;
      t_trxframe _fr;
  
      snew(_oenv, 1);
      output_env_init_default(_oenv);
  
      if(not read_first_frame(_oenv, &_status, trjName.c_str(), &_fr, 0))
      {
        output_env_done(_oenv);
        return 0;
      }  
      nFrames = 1;

      while(read_next_frame(_oenv, _status, &_fr))
        nFrames++;

      close_trx(_status);
      output_env_done(_oenv);
    }
    else
    	nFrames = getLastFrameTime() / getTimeStep() + 1;

    return cachedNFrames = nFrames;
  }

  unsigned int
  Gromacs::getCurrentFrame() const
  {
    return currentFrame + 1;
  }

  void
  Gromacs::waitNextFrame() const
  {
    operationMutex.lock();
    if(not abortFlag and operationThread.joinable()
       and getCurrentFrame() < getFramesCount())
    {
      std::unique_lock<std::mutex> slock(operationMutex, std::defer_lock);
      wakeCondition.wait(slock);
    }
    operationMutex.unlock();
  }

  void
  Gromacs::waitNextFrame(unsigned int refFrame) const
  {
    operationMutex.lock();
    while(not abortFlag and operationThread.joinable()
          and getCurrentFrame() < refFrame + 1)
    {
      std::unique_lock<std::mutex> slock(operationMutex, std::defer_lock);
      wakeCondition.wait(slock);
    }
    operationMutex.unlock();
  }

  float
  Gromacs::getLastFrameTime() const
  {
    if(not exists(trjName))
      return 0;

    output_env_t _oenv;
    t_trxstatus *_status;
    t_trxframe _fr;
    t_fileio *fio;
    gmx_bool ok;

    snew(_oenv, 1);
    output_env_init_default(_oenv);

    if(not read_first_frame(_oenv, &_status, trjName.c_str(), &_fr, 0))
    {
      output_env_done(_oenv);
      return 0;
    }

    fio = trx_get_fileio(_status);
    float lastTime = xdr_xtc_get_last_frame_time(gmx_fio_getfp(fio),
                                                 gmx_fio_getxdr(fio),
                                                 _fr.natoms, &ok);

    close_trx(_status);
    output_env_done(_oenv);
    if(not ok)
      throw "ERROR: can not get the last frame time";

    return lastTime;
  }

  float
  Gromacs::getTimeStep() const
  {
    if(not exists(trjName))
      return 0;

    if(timeStepCached != 0)
      return timeStepCached;

    output_env_t _oenv;
    t_trxstatus *_status;
    t_trxframe _fr;
    float time1, time2;

    snew(_oenv, 1);
    output_env_init_default(_oenv);

    if(not read_first_frame(_oenv, &_status, trjName.c_str(), &_fr, 0))
    {
      output_env_done(_oenv);
      return 0;
    }

    time1 = _fr.time;
    read_next_frame(_oenv, _status, &_fr);
    time2 = _fr.time;
    close_trx(_status);

    output_env_done(_oenv);

    timeStepCached = time2 - time1;
    return timeStepCached;
  }

  unsigned int
  Gromacs::getFrameStep() const
  {
    if(exists(trjName))
      return 0;

    output_env_t _oenv;
    t_trxstatus *_status;
    t_trxframe _fr;

    snew(_oenv, 1);
    output_env_init_default(_oenv);

    if(read_first_frame(_oenv, &_status, trjName.c_str(), &_fr, 0))
    {
      close_trx(_status);
      output_env_done(_oenv);
      return 0;
    }

    read_next_frame(_oenv, _status, &_fr);
    close_trx(_status);

    output_env_done(_oenv);

    return _fr.step;
  }

  float
  Gromacs::getBegin() const
  {
    return _begin;
  }

  void
  Gromacs::setBegin(float beginTime)
  {
    _begin = beginTime;
    setTimeValue(TBEGIN, beginTime);
  }

  float
  Gromacs::getEnd() const
  {
    return _end;
  }

  void
  Gromacs::setEnd(float endTime)
  {
    _end = endTime;
    setTimeValue(TEND, endTime);
  }

  void
  Gromacs::abort()
  {
    abortFlag = true;
    if(operationThread.joinable())
      operationThread.join();
    wakeCondition.notify_all();
  }

  bool
  Gromacs::isAborting() const
  {
    return abortFlag;
  }

  bool
  Gromacs::usePBC() const noexcept
  {
    return _usePBC;
  }

  bool
  Gromacs::usePBC(bool value) noexcept
  {
    return _usePBC = value;
  }

  template void Gromacs::__calculateSas(Session<std::fstream>&);
  template void Gromacs::__calculateSas(Session<std::ofstream>&);

  template void Gromacs::calculateSas(Session<std::fstream>&);
  template void Gromacs::calculateSas(Session<std::ofstream>&);
}
