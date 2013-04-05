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

#ifndef _SASANALYSIS_H
#define _SASANALYSIS_H

namespace PstpFinder
{
  template<typename T> class SasAnalysis_Base;
  template<typename T> class SasAnalysis_Read;
  template<typename T> class SasAnalysis_Write;
  template<typename T, typename = void> class SasAnalysis;
}

#include "utils.h"
#include "SasAnalysisThread.h"
#include "Serializer.h"
#include "Gromacs.h"
#include "Session.h"
#include "SasAtom.h"

#include <thread>
#include <vector>
#include <boost/circular_buffer.hpp>

namespace PstpFinder
{
  template<typename T, typename>
  class SasAnalysisThread;

  template<typename T>
  class SasAnalysis_Base
  {
    public:
      SasAnalysis_Base(unsigned int nAtoms, const Gromacs&, Session<T>&);
      SasAnalysis_Base(const Gromacs& gromacs, const string& sessionFileName);
      SasAnalysis_Base(const Gromacs&, Session<T>&);
      virtual bool setMaxBytes(unsigned long bytes);
      virtual unsigned long getMaxBytes();
      virtual bool setMaxChunkSize(unsigned long bytes);
      virtual unsigned long getMaxChunkSize();

    protected:
      boost::circular_buffer<std::vector<vector<SasAtom>>> chunks;
      vector<vector<SasAtom>> frames;
      unsigned int nAtoms;
      Session<T> rawSession;
      MetaStream<T>& sasMetaStream;
      Serializer<MetaStream<T>>* serializer;
      std::streampos fileStreamEnd;
      const Gromacs* gromacs;
      unsigned long maxFrames, maxBytes, maxChunk;
      unsigned int bufferCount;
      unsigned int bufferMax;
      mutable condition_variable bufferCountCondition;
      mutable mutex bufferMutex, bufferCountMutex;
      bool changeable;

      template<typename, typename> friend class SasAnalysisThread_Base;
      template<typename, typename> friend class SasAnalysisThread;
      typedef SasAnalysisThread<T> SasAnalysisThreadType;
      SasAnalysisThreadType* analysisThread;

      virtual void init();
      virtual void updateChunks();
  };

  template<typename T>
  class SasAnalysis_Read : public SasAnalysis_Base<T>
  {
    public:
      typedef SasAnalysisThread<T> SasAnalysisThreadType;
      SasAnalysis_Read(unsigned int nAtoms, const Gromacs& gromacs,
                        Session<T>& sessionFile) :
          Base(nAtoms, gromacs, sessionFile),
          beginTriggered(false){ updateChunks(); }
      SasAnalysis_Read(const Gromacs& gromacs, const string& sessionFileName) :
          Base(gromacs, sessionFileName),
          beginTriggered(false) { updateChunks(); }
      SasAnalysis_Read(const Gromacs& gromacs, Session<T>& sessionFile) :
          Base(gromacs, sessionFile),
          beginTriggered(false){ updateChunks(); }
      virtual ~SasAnalysis_Read();
      virtual SasAnalysis_Read& removeMe(const vector<SasAtom>*& sasAtom);

      class const_iterator :
          public std::iterator<std::forward_iterator_tag, const vector<SasAtom>>
      {
          friend class SasAnalysis_Read<T>;

        public:
          reference operator*() const;
          pointer operator->() const;
          const_iterator& operator++();
          bool operator!=(const const_iterator& other) const;
          const_iterator() : parent(nullptr), iter() {}

        private:
          const_iterator(SasAnalysis_Read<T>* parent)
            : parent(parent), iter() { handleEmptyFrames(); }
          const_iterator(SasAnalysis_Read<T>* parent,
                         vector<vector<SasAtom>>::const_iterator iter)
            : parent(parent), iter(iter) {}
          void handleEmptyFrames();

          SasAnalysis_Read<T>* parent;
          vector<vector<SasAtom>>::const_iterator iter;
      };
      friend class const_iterator;

      const_iterator begin();
      const_iterator end();

    private:
      typedef SasAnalysis_Base<T> Base;
      template<typename, typename> friend class SasAnalysisThread_Base;
      template<typename, typename> friend class SasAnalysisThread;

      virtual vector<vector<SasAtom>>
        loadChunk(Serializer<MetaStream<T>>& in);
      bool open();
      virtual void updateChunks();

      bool beginTriggered;
  };


  template<typename T>
  class SasAnalysis_Write : public SasAnalysis_Base<T>
  {
    public:
      typedef SasAnalysisThread<T> SasAnalysisThreadType;
      SasAnalysis_Write(unsigned int nAtoms, const Gromacs& gromacs,
                        Session<T>& sessionFile) :
          Base(nAtoms, gromacs, sessionFile), readFrames(0) { updateChunks(); }
      SasAnalysis_Write(const Gromacs& gromacs, const string& sessionFileName) :
          Base(gromacs, sessionFileName), readFrames(0) { updateChunks(); }
      SasAnalysis_Write(const Gromacs& gromacs, Session<T>& sessionFile) :
          Base(gromacs, sessionFile), readFrames(0) { updateChunks(); }
      virtual ~SasAnalysis_Write();
      virtual const SasAnalysis_Write& push_back(const vector<SasAtom>& sasAtoms);
      unsigned int getReadFrames() const;

    protected:
      unsigned int readFrames;

    private:
      typedef SasAnalysis_Base<T> Base;
      template<typename, typename> friend class SasAnalysisThread_Base;
      template<typename, typename> friend class SasAnalysisThread;

      virtual void dumpChunk(const vector<vector<SasAtom>>& chunk,
                Serializer<MetaStream<T>>& out) const;
      virtual void flush();
      virtual bool save();
      virtual void updateChunks();
  };

  template<typename T>
  class SasAnalysis<T,
      typename enable_if<is_base_of<base_stream(basic_istream, T), T>::value
                          and not is_base_of<base_stream(basic_ostream, T),
                                              T>::value>
                ::type> :
      public SasAnalysis_Read<T>
  {
    public:
      SasAnalysis(unsigned int nAtoms, const Gromacs& gromacs,
                  Session<T>& sessionFile) :
          SasAnalysis_Read<T>(nAtoms, gromacs, sessionFile) {}
      SasAnalysis(const Gromacs& gromacs, Session<T>& sessionFile) :
          SasAnalysis_Read<T>(gromacs, sessionFile) {}
      SasAnalysis(const Gromacs& gromacs, const string& sessionFileName) :
          SasAnalysis_Read<T>(gromacs, sessionFileName) {}
  };
}

#include "SasAnalysis.cpp"
#endif /* _SASANALYSIS_H */
