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

#include "Pdb.h"

#include <string>
#include <thread>
#include <mutex>
#include <condition_variable>

#if GMXVER <= 45
/* Workaround - is not defined as "C", let's include it before others */
extern "C"
{
#include <gromacs/gmxfio.h>
}

#include <gromacs/atomprop.h>
#include <gromacs/statutil.h>

#else
#include <gromacs/legacyheaders/atomprop.h>
#include <gromacs/fileio/trxio.h>
#endif /* GMXVER <= 45 */

namespace gmx_legacy
{
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

#if GMXVER == 50

/* Taken from src/gromacs/fileio/gmx_system_xdr.h, Gromacs 5.0.1 */
typedef int bool_t;
typedef int xdr_int32_t;
typedef unsigned int xdr_uint32_t;

enum xdr_op {
    XDR_ENCODE = 0,
    XDR_DECODE = 1,
    XDR_FREE   = 2
};

struct XDR
{
    enum xdr_op x_op;       /* operation; fast additional param */
    struct xdr_ops
    {
        bool_t       (*x_getbytes) (XDR *__xdrs, char *__addr, unsigned int __len);
        /* get some bytes from " */
        bool_t       (*x_putbytes) (XDR *__xdrs, char *__addr, unsigned int __len);
        /* put some bytes to " */
        unsigned int (*x_getpostn) (XDR *__xdrs);
        /* returns bytes off from beginning */
        bool_t       (*x_setpostn) (XDR *__xdrs, unsigned int __pos);
        /* lets you reposition the stream */
        xdr_int32_t *(*x_inline) (XDR *__xdrs, int __len);
        /* buf quick ptr to buffered data */
        void         (*x_destroy) (XDR *__xdrs);
        /* free privates of this xdr_stream */
        bool_t       (*x_getint32) (XDR *__xdrs, xdr_int32_t *__ip);
        /* get a int from underlying stream */
        bool_t       (*x_putint32) (XDR *__xdrs, xdr_int32_t *__ip);
        /* put a int to " */
        bool_t       (*x_getuint32) (XDR *__xdrs, xdr_uint32_t *__ip);
        /* get a unsigned int from underlying stream */
        bool_t       (*x_putuint32) (XDR *__xdrs, xdr_uint32_t *__ip);
        /* put a int to " */
    }
    *x_ops;
    char *x_public;     /* users' data */
    char *x_private;    /* pointer to private data */
    char *x_base;       /* private used for position info */
    int   x_handy;      /* extra private word */
};

#endif /* GMXVER == 50 */

/* Taken from src/gromacs/fileio/xdrf.h, Gromacs 5.0.1 */
extern "C" float xdr_xtc_get_last_frame_time(FILE *fp, XDR *xdrs, int natoms, gmx_bool * bOK);
extern "C" XDR *gmx_fio_getxdr(struct t_fileio *fio);

/* Taken from src/gromacs/fileio/timecontrol.h, Gromacs 5.0.1 */
enum {
    TBEGIN, TEND, TDELTA, TNR
};

extern "C" void setTimeValue(int tcontrol, real value);

} /* namespace gmx_legacy */

namespace PstpFinder
{
  class Gromacs
  {
    public:
      Gromacs(float solventSize = 0.14);
      Gromacs(const std::string& trajectoryFileName,
              const std::string& topologyFileName,
              float solventSize = 0.14);
      Gromacs(const Gromacs& gromacs);
      ~Gromacs();
      // FIXME: It will have to return an object Molecular Dynamics with SAS
      // FIXME: additional informations
      template<typename Stream>
      void calculateSas(Session<Stream>& session);
    
      std::string getTrajectoryFile() const;
      std::string getTopologyFile() const;
      unsigned long getAtomsCount() const;
      unsigned int getFramesCount() const;
      float getTimeStep() const; // In nsec
      float getLastFrameTime() const;
      unsigned int getFrameStep() const;
      float getBegin() const;
      void setBegin(float beginTime);
      float getEnd() const;
      void setEnd(float endTime);
      std::vector<atom_id> getGroup(const std::string& groupName);
      std::vector<atom_id> getGroup(const std::string& groupName) const;
    
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
      const Protein<>& __calculateAverageStructure();
      void calculateAverageStructure();
      const Protein<>& getAverageStructure() const;
      void setAverageStructure(Protein<> structure);
      void waitOperation();
      void abort();
      bool isAborting() const;
      bool usePBC() const noexcept;
      bool usePBC(bool value) noexcept;

    private:
#if GMXVER >= 45
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
      gmx_mtop_t mtop;
      bool _usePBC;

      std::string trjName, tprName;
      bool gotTrajectory, gotTopology, readyToGetX;
      float solSize;
      std::string sasTarget;
    
      std::thread operationThread;
      mutable std::mutex operationMutex;
      mutable std::condition_variable wakeCondition;
      mutable unsigned int cachedNFrames;
      unsigned int currentFrame; // index-0 based -- like always
      Protein<> averageStructure;
      float _begin, _end;
      mutable float timeStepCached;
      bool abortFlag;
    
      void init(float solventSize);
      bool getTopology();
      bool getTrajectory();
      bool readNextX();
  };
}
#endif

