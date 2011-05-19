#ifndef _GROMACS_H
#define _GROMACS_H

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "Protein.h"

#include <string>

#include <boost/thread/thread.hpp>
#include <boost/interprocess/sync/interprocess_mutex.hpp>
#include <boost/interprocess/sync/interprocess_condition.hpp>

#include <atomprop.h>
#include <statutil.h>

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
    void calculateSas();
    
    string getTrajectoryFile() const;
    string getTopologyFile() const;
    unsigned long getAtomsCount() const;
    unsigned int getFramesCount() const;
    float getTimeStep() const; // In nsec
    unsigned int getFrameStep() const;
    void setBegin(float beginTime);
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

    void __calculateSas();
    const Protein& __calculateAverageStructure();
    const Protein& getAverageStructure() const;
    // TODO: calculateAverageStructureDetach()
    void waitOperation();
  private:
#ifdef GMX45
    output_env_t oenv;
    t_trxstatus *status;
#else
    t_commrec *cr;
    int status, step;
    real lambda;
#endif
    real t;
    rvec *x, *xtop;
    matrix box;
    int natoms, ePBC;
    t_topology top;

    gmx_atomprop_t aps;
    string trjName, tprName;
    bool gotTrajectory, gotTopology;
    float solSize;
    string sasTarget;
    gmx_mtop_t mtop;
    
    boost::thread operationThread;
    mutable boost::interprocess::interprocess_mutex operationMutex, wakeMutex;
    mutable boost::interprocess::interprocess_condition wakeCondition;
    mutable unsigned int cachedNFrames;
    unsigned int currentFrame; // index-0 based -- like always
    Protein averageStructure;
    
    void init(float solventSize);
    bool getTopology();
    bool getTrajectory();
    bool readNextX();
  };
};

#endif

