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

#ifndef _OPERATIONTHREAD_H
#define _OPERATIONTHREAD_H

#include "SasAnalysis.h"

namespace PstpFinder
{
  template<typename T>
  class SasAnalysis;

  template<typename T>
  class OperationThread_Base
  {
    public:
      OperationThread_Base(SasAnalysis<T>& parent);
      virtual ~OperationThread_Base();
      virtual void wakeUp();
      virtual void stop();
      void threadSave();
      void threadOpen();

    protected:
      SasAnalysis<T>* parent;
      bool isStopped;
      thread operationThread;
      condition_variable wakeCondition;
      mutex wakeMutex;
  };

  template<typename T>
  OperationThread_Base<T>::OperationThread_Base(SasAnalysis<T>& parent)
  {
    this->parent = &parent;
    isStopped = false;
  }

  template<typename T>
  OperationThread_Base<T>::~OperationThread_Base()
  {
    stop();
  }

  template<typename T>
  void
  OperationThread_Base<T>::threadSave()
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

  template<typename T>
  void
  OperationThread_Base<T>::threadOpen()
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

  template<typename T>
  void
  OperationThread_Base<T>::wakeUp()
  {
    wakeCondition.notify_one();
  }

  template<typename T>
  void
  OperationThread_Base<T>::stop()
  {
    if(isStopped)
    {
      if(operationThread.joinable())
        operationThread.join();
      return;
    }

    parent->bufferMutex.lock();
    isStopped = true;
    wakeCondition.notify_one();
    parent->bufferMutex.unlock();

    if(operationThread.joinable())
      operationThread.join();
  }

  template<typename T, typename = void>
  class OperationThread;

  template<typename T>
  class OperationThread<T, typename enable_if<
    is_base_of<base_stream(basic_istream), T>::value and
    not is_base_of<base_stream(basic_ostream), T>::value>::type>
    : public OperationThread_Base<T>
  {
    public:
      OperationThread(SasAnalysis<T>& parent) :
        OperationThread_Base<T>(parent)
      {
        this->operationThread = thread(
              bind(&OperationThread_Base<T>::threadOpen, ref(*this)));
      }
      virtual ~OperationThread() {}
      virtual void threadOpen() { OperationThread_Base<T>::threadOpen(); }
  };

  template<typename T>
  class OperationThread<T, typename enable_if<
    not is_base_of<base_stream(basic_istream), T>::value and
    is_base_of<base_stream(basic_ostream), T>::value>::type>
    : public OperationThread_Base<T>
  {
    public:
      OperationThread(SasAnalysis<T>& parent):
        OperationThread_Base<T>(parent)
      {
        this->operationThread = thread(
            bind(&OperationThread_Base<T>::threadSave, ref(*this)));
      }
      virtual ~OperationThread() {}
      virtual void threadSave() { OperationThread_Base<T>::threadSave(); }
  };

  template<typename T>
  class OperationThread<T, typename enable_if<
    is_base_of<base_stream(basic_istream), T>::value and
    is_base_of<base_stream(basic_ostream), T>::value>::type>
    : public OperationThread_Base<T>
  {
    public:
      OperationThread(SasAnalysis<T>& parent):
        OperationThread_Base<T>(parent)
      {
        if(parent.mode == SasAnalysis<T>::MODE_OPEN)
        {
          this->operationThread = thread(
                bind(&OperationThread_Base<T>::threadOpen, ref(*this)));
        }
        else
        {
          this->operationThread = thread(
              bind(&OperationThread_Base<T>::threadSave, ref(*this)));
        }
      }
      virtual ~OperationThread() {}
      virtual void threadSave() { OperationThread_Base<T>::threadSave(); }
  };
}

#endif /* _OPERATIONTHREAD_H */
