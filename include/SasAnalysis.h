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

#include "SasAtom.h"
#include "Gromacs.h"
#include "Session.h"

#include <boost/thread.hpp>
#include <boost/circular_buffer.hpp>
#include <boost/filesystem.hpp>
#include <boost/thread/detail/thread.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/device/file_descriptor.hpp>
#include <boost/iostreams/filter/zlib.hpp>
#include <boost/range/iterator_range.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/interprocess/sync/interprocess_semaphore.hpp>
#include <boost/interprocess/sync/interprocess_mutex.hpp>
#include <boost/interprocess/sync/interprocess_condition.hpp>
#include <boost/interprocess/sync/scoped_lock.hpp>

#include <vector>
#include <string>
#include <iostream>
#include <sstream>

namespace PstpFinder
{
  class SasAnalysis
  {
    public:
      SasAnalysis(unsigned int nAtoms, std::string filename = "/tmp/sas.psf",
                  bool savingMode = true);

      SasAnalysis(const Gromacs& gromacs, std::string filename = "/tmp/sas.psf",
                  bool savingMode = true);

      ~SasAnalysis();

      const SasAnalysis&
      operator <<(SasAtom* sasAtoms);

      SasAnalysis&
      operator >>(SasAtom*& sasAtom);

      bool
      setMaxBytes(unsigned long bytes);

      unsigned long
      getMaxBytes();

      bool
      setMaxChunkSize(unsigned long bytes);

      unsigned long
      getMaxChunkSize();

    private:
      boost::circular_buffer<std::vector<SasAtom*> > chunks;
      std::vector<SasAtom*> frames;
      unsigned int nAtoms;
      std::string filename;
      boost::iostreams::file_descriptor fileIO;
      Session sessionFile;
      std::streampos fileStreamEnd;
      const Gromacs* gromacs;
      unsigned long maxFrames, maxBytes, maxChunk;
      mutable boost::interprocess::interprocess_semaphore* bufferSemaphore;
      unsigned int bufferSemaphoreCount;
      unsigned int bufferSemaphoreMax;
      mutable boost::interprocess::interprocess_mutex bufferMutex;
      bool changeable;
      boost::iostreams::filtering_istream inFilter;
      boost::archive::binary_iarchive* inArchive;
      boost::iostreams::filtering_ostream outFilter;
      boost::archive::binary_oarchive* outArchive;

      enum
      {
        MODE_OPEN,
        MODE_SAVE
      } mode;

      class OperationThread
      {
        public:
          OperationThread(SasAnalysis& parent);
          ~OperationThread();

          void
          wakeUp();

          void
          stop();

          void
          threadSave();

          void
          threadOpen();

        private:
          SasAnalysis* parent;
          bool isStopped;
          boost::thread thread;
          boost::interprocess::interprocess_condition wakeCondition;
          boost::interprocess::interprocess_mutex wakeMutex;
      };

      friend class OperationThread;
      OperationThread* operationThread;

      void
      init(string saveFile, bool savingMode = true);

      void
      dumpChunk(const std::vector<SasAtom*>& chunk,
                boost::archive::binary_oarchive& out) const;

      std::vector<SasAtom*>
      loadChunk(boost::archive::binary_iarchive& in);

      void
      flush();

      bool
      save(const std::string& filename);

      bool
      save();

      bool
      open();

      void
      updateChunks();
  };
}
;

#endif

