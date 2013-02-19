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

#ifndef _SASANALYSISTHREAD_H
#define _SASANALYSISTHREAD_H

#include "SasAnalysis.h"

#include <thread>

using namespace std;

namespace PstpFinder
{
  template<typename T, typename SasClass>
  class SasAnalysisThread_Base
  {
    public:
      SasAnalysisThread_Base(SasClass& parent);
      virtual ~SasAnalysisThread_Base();
      virtual void wakeUp();
      virtual void stop();
      void threadSave();
      void threadOpen();

    protected:
      SasClass* parent;
      bool isStopped;
      thread analysisThread;
      condition_variable wakeCondition;
      mutex wakeMutex;
  };

  template<typename T, typename SasClass>
  SasAnalysisThread_Base<T, SasClass>::SasAnalysisThread_Base(SasClass& parent)
  {
    this->parent = &parent;
    isStopped = false;
  }

  template<typename T, typename SasClass>
  SasAnalysisThread_Base<T, SasClass>::~SasAnalysisThread_Base()
  {
    stop();
  }

  template<typename T, typename SasClass>
  void
  SasAnalysisThread_Base<T, SasClass>::threadSave()
  {
    while(not isStopped)
    {
      parent->bufferMutex.lock();

      while(parent->bufferCount == parent->bufferMax - 1 and not isStopped)
      {
        unique_lock<mutex> lock(wakeMutex);
        parent->bufferMutex.unlock();
        wakeCondition.wait(lock);
        parent->bufferMutex.lock();
      }

      if(isStopped or (parent->gromacs and parent->gromacs->isAborting()))
      {
        parent->bufferMutex.unlock();
        break;
      }

      if(parent->save())
      {
        vector<SasAtom*>& curChunk = parent->chunks.front();
        for(std::vector<SasAtom*>::iterator i = curChunk.begin();
            i < curChunk.end(); i++)
          delete[] *i;
        curChunk.clear();

        parent->chunks.pop_front();
        parent->bufferCount++;
        parent->bufferCountCondition.notify_all();
      }

      parent->bufferMutex.unlock();
    }

    if(parent->gromacs and parent->gromacs->isAborting())
      return;

    parent->bufferMutex.lock();
    while(parent->bufferCount < parent->bufferMax - 1)
    {
      parent->save();
      parent->chunks.pop_front();
      parent->bufferCount++;
      parent->bufferCountCondition.notify_all();
    }
    parent->bufferMutex.unlock();
  }

  template<typename T, typename SasClass>
  void
  SasAnalysisThread_Base<T, SasClass>::threadOpen()
  {
    while(not isStopped)
    {
      parent->bufferMutex.lock();

      while(parent->bufferCount == parent->bufferMax - 1 and not isStopped)
      {
        unique_lock<mutex> lock(wakeMutex);
        parent->bufferMutex.unlock();
        wakeCondition.wait(lock);
        parent->bufferMutex.lock();
      }

      if(isStopped)
      {
        while(parent->bufferCount < parent->bufferMax - 1)
        {
          parent->bufferCount++;
          parent->bufferCountCondition.notify_all();
        }
        parent->bufferMutex.unlock();
        break;
      }

      if(parent->gromacs and parent->gromacs->isAborting())
      {
        parent->bufferMutex.unlock();
        break;
      }

      if(not parent->open())
        isStopped = true;

      parent->bufferCount++;
      parent->bufferCountCondition.notify_all();
      parent->bufferMutex.unlock();
    }
  }

  template<typename T, typename SasClass>
  void
  SasAnalysisThread_Base<T, SasClass>::wakeUp()
  {
    wakeCondition.notify_one();
  }

  template<typename T, typename SasClass>
  void
  SasAnalysisThread_Base<T, SasClass>::stop()
  {
    if(isStopped)
    {
      if(analysisThread.joinable())
        analysisThread.join();
      return;
    }

    parent->bufferMutex.lock();
    isStopped = true;
    wakeCondition.notify_one();
    parent->bufferMutex.unlock();

    if(analysisThread.joinable())
      analysisThread.join();
  }

  template<typename T, typename = void>
  class SasAnalysisThread;

  template<typename T>
  class SasAnalysisThread<T, typename enable_if<
    is_base_of<base_stream(basic_istream, T), T>::value and
    not is_base_of<base_stream(basic_ostream, T), T>::value>::type>
    : public SasAnalysisThread_Base<T, SasAnalysis_Read<T>>
  {
    public:
      typedef SasAnalysisThread_Base<T, SasAnalysis_Read<T>> Base;
      SasAnalysisThread(SasAnalysis_Read<T>& parent) :
        Base(parent)
      {
        this->analysisThread = thread(
              bind(&Base::threadOpen, ref(*this)));
      }
      virtual ~SasAnalysisThread() {}
      virtual void threadOpen() { Base::threadOpen(); }
  };

  template<typename T>
  class SasAnalysisThread<T, typename enable_if<
    is_base_of<base_stream(basic_ostream, T), T>::value>::type>
    : public SasAnalysisThread_Base<T, SasAnalysis_Write<T>>
  {
    public:
      typedef SasAnalysisThread_Base<T, SasAnalysis_Write<T>> Base;
      SasAnalysisThread(SasAnalysis_Write<T>& parent):
        Base(parent)
      {
        this->analysisThread = thread(
            bind(&Base::threadSave, ref(*this)));
      }
      virtual ~SasAnalysisThread() {}
      virtual void
      threadSave() { Base::threadSave(); }

      void waitForFlush()
      {
        if(Base::isStopped)
          return;

        Base::parent->bufferMutex.lock();
        while(Base::parent->bufferCount < Base::parent->bufferMax - 1 and not Base::isStopped)
        {
          unique_lock<mutex> lock(Base::wakeMutex);
          Base::parent->bufferMutex.unlock();
          Base::wakeCondition.wait(lock);
          Base::parent->bufferMutex.lock();

          if(Base::parent->gromacs and Base::parent->gromacs->isAborting())
          {
            Base::parent->bufferMutex.unlock();
            break;
          }
        }
        Base::parent->bufferMutex.unlock();
      }
  };
}

#endif /* _SASANALYSISTHREAD_H */
