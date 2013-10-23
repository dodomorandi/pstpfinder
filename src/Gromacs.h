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

// This must be included BEFORE defining _H scope
#include "Session.h"

#ifndef _GROMACS_H
#define _GROMACS_H

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "Protein.h"

#include <string>
#include <thread>
#include <mutex>
#include <condition_variable>

/* Workaround - is not defined as "C", let's include it before others */
extern "C"
{
#include <gromacs/gmxfio.h>
}

#include <gromacs/atomprop.h>
#include <gromacs/statutil.h>

/* Taken from src/tools/nsc.h, Gromacs 4.5.2 */
#define FLAG_ATOM_AREA  04
extern "C" int
    nsc_dclm_pbc(rvec *coords, real *radius, int nat, int densit, int mode,
                 real *value_of_area, real **at_area, real *value_of_vol,
                 real **lidots, int *nu_dots, atom_id index[], int ePBC,
                 matrix box);

/* Taken from src/gmxlib/gmx_atomprop.c, Gromacs 4.5.2 */
typedef struct {
  gmx_bool   bSet;
  int    nprop,maxprop;
  char   *db;
  double def;
  char   **atomnm;
  char   **resnm;
  gmx_bool   *bAvail;
  real   *value;
} aprop_t;

struct gmx_atomprop {
  gmx_bool       bWarned;
  aprop_t    prop[epropNR];
  gmx_residuetype_t restype;
};

/* Taken from src/gmxlib/index.c, Gromacs 4.5.2 */
struct gmx_residuetype
{
    int      n;
    char **  resname;
    char **  restype;
};

using namespace std;

namespace PstpFinder
{
  class Gromacs
  {
    public:
      Gromacs(float solventSize = 0.14);
      Gromacs(const string& trajectoryFileName, const string& topologyFileName,
              float solventSize = 0.14);
      Gromacs(const Gromacs& gromacs);
      ~Gromacs();
      // FIXME: It will have to return an object Molecular Dynamics with SAS
      // FIXME: additional informations
      template<typename Stream>
      void calculateSas(Session<Stream>& session);
    
      string getTrajectoryFile() const;
      string getTopologyFile() const;
      unsigned long getAtomsCount() const;
      unsigned int getFramesCount() const;
      float getTimeStep() const; // In nsec
      unsigned int getFrameStep() const;
      float getBegin() const;
      void setBegin(float beginTime);
      float getEnd() const;
      void setEnd(float endTime);
      vector<atom_id> getGroup(const string& groupName);
      vector<atom_id> getGroup(const string& groupName) const;
    
      /**
       * @brief Returns the number of the current frame starting from 1
       */
      unsigned int getCurrentFrame() const;
      void waitNextFrame() const;

      /**
       * @brief Waits for next frame
       * @param refFrame Index-1 based of the hypotetical current frame.
       *                 The method will return when the current frame is
       *                 refFrame + 1
       */
      void waitNextFrame(unsigned int refFrame) const;

      template<typename Stream>
      void __calculateSas(Session<Stream>& session);
      const Protein& __calculateAverageStructure();
      void calculateAverageStructure();
      const Protein& getAverageStructure() const;
      void setAverageStructure(Protein structure);
      void waitOperation();
      void abort();
      bool isAborting() const;
      bool usePBC() const noexcept;
      bool usePBC(bool value) noexcept;
    private:
#ifdef GMX45
      output_env_t oenv;
      t_trxstatus *status;
#else
      t_commrec *cr;
      int status, step;
      real lambda;
#endif
      rvec *xtop;
      int natoms, ePBC;
      t_topology top;
      t_trxframe fr;

      gmx_atomprop_t aps;
      string trjName, tprName;
      bool gotTrajectory, gotTopology, readyToGetX;
      float solSize;
      string sasTarget;
      gmx_mtop_t mtop;
    
      thread operationThread;
      mutable mutex operationMutex;
      mutable condition_variable wakeCondition;
      mutable unsigned int cachedNFrames;
      unsigned int currentFrame; // index-0 based -- like always
      Protein averageStructure;
      float _begin, _end;
      mutable float timeStepCached;
      bool abortFlag;
      bool _usePBC;
    
      void init(float solventSize);
      bool getTopology();
      bool getTrajectory();
      bool readNextX();
  };
}
#endif

