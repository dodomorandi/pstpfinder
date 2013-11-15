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

#ifndef _PITTPI_H
#define _PITTPI_H

#include "config.h"
#include "Gromacs.h"
#include "Protein.h"
#include "SasAtom.h"

#include <vector>
#include <thread>

namespace PstpFinder
{
  class Group
  {
    public:
      Group(const Residue<SasPdbAtom>& refResidue);
      Group(const SasPdbAtom& refAtomH);
      Group(const SasPdbAtom& refAtomH, const Protein<SasPdbAtom>& protein);
      Group(const Group& group);
      Group(Group&& group);
      Group(const Group& group, const Protein<SasPdbAtom>& protein);
      Group& operator =(const Group& group);
      Group& operator =(Group&& group);
      Group& operator <<(const Residue<SasPdbAtom>& value);
      Group& operator <<(const Group& group);
      const std::vector<const Residue<SasPdbAtom>*>& getResidues() const;
      const SasPdbAtom& getCentralH() const;
      const Residue<SasPdbAtom>& getCentralRes() const;

      std::vector<float> sas;
      unsigned int zeros;
      protected:
      Group() : referenceAtom(nullptr), referenceRes(nullptr)
      {}

      const SasPdbAtom* referenceAtom;
      const Residue<SasPdbAtom>* referenceRes;
      std::vector<const Residue<SasPdbAtom>*> residues;
    };

  struct Pocket
  {
      const Group* group;
      unsigned int startFrame;
      float startPs;
      unsigned int endFrame;
      float endPs;
      unsigned int width;
      float openingFraction;
      unsigned int averageNearFrame;
      float averageNearPs;
      unsigned int maxAreaFrame;
      float maxAreaPs;

      Pocket(const Group& group) :
          group(&group)
      {
      }
      Pocket(const Pocket& pocket) :
          group(pocket.group)
      {
        *this << pocket;
      }

      Pocket&
      operator =(const Pocket& pocket)
      {
        group = pocket.group;
        *this << pocket;
        return *this;
      }

      Pocket&
      operator <<(const Pocket& pocket)
      {
        startFrame = pocket.startFrame;
        startPs = pocket.startPs;
        endFrame = pocket.endFrame;
        endPs = pocket.endPs;
        width = pocket.width;
        openingFraction = pocket.openingFraction;
        averageNearFrame = pocket.averageNearFrame;
        averageNearPs = pocket.averageNearPs;
        maxAreaFrame = pocket.maxAreaFrame;
        maxAreaPs = pocket.maxAreaPs;

        return *this;
      }

    protected:
      Pocket() :
          group(nullptr)
      {
      }
  };

  /**
   * @brief Related to the method developed by Matteo De Chiara and Silvia
   * Bottini
   *
   * The original code have been developed in Perl. A code refactory have been
   * done to obtain a more usable and stable program.
   *
   * @note This is only a stub. The class must be rewritten, but for now
   * the main objective remains the productivity of the software. I created
   * a class only to remember what must be done in the future. Possibly soon.
   * @author Edoardo Morandi
   */
  class Pittpi
  {
    public:
      Pittpi(Gromacs& gromacs, const std::string& sessionFileName, float radius,
             unsigned long threshold, bool runPittpi = true);
      Pittpi(const Pittpi& pittpi);
      Pittpi(const Pittpi& pittpi, const Gromacs& gromacs);
      ~Pittpi();

      void
      join();
      void
      setStatus(float value) const;
      void
      setStatusDescription(const std::string& description) const;
      float
      getStatus() const;
      std::string
      getStatusDescription() const;
      void
      waitNextStatus();
      bool
      isFinished();
      void
      abort();
      const std::vector<Pocket>&
      getPockets() const;

      template<typename Stream>
        void
        save(Stream& stream) const;
      template<typename Stream>
        void
        load(Stream& stream);
    private:
      class SerializablePockets
      {
        public:
          struct SerializablePocket : public Pocket
          {
              template<typename, typename >
                friend class Serializer;
              long groupIndex;

              SerializablePocket() = default;
              SerializablePocket(const Pocket& pocket) : Pocket(pocket)
              {}

              template<typename Serializer>
              void serialize(Serializer serializer);
            };

            SerializablePockets() = default;
            SerializablePockets(const std::vector<Pocket>& pockets,
                const std::vector<Group>& groups);

            template<typename Serializer>
            void serialize(Serializer serializer);
            void updatePockets(std::vector<Pocket>& pocketsToUpdate,
                const std::vector<Group>& groups) const;

            private:
            std::vector<SerializablePocket> pockets;
          };

          class SerializableGroups
          {
            public:
            class SerializableGroup : public Group
            {
              public:
              template<typename, typename> friend class Serializer;
              friend class Pittpi::SerializableGroups;
              long referenceAtomIndex, referenceResIndex;
              std::vector<long> residuesIndex;

              SerializableGroup() = default;
              SerializableGroup(const Group& group) : Group(group)
              {}

              template<typename Serializer>
              void serialize(Serializer serializer);
            };

            SerializableGroups() = default;
            SerializableGroups(const std::vector<Group>& groups,
                const Protein<SasPdbAtom>& protein);

            template<typename Serializer>
            void serialize(Serializer serializer);
            void updateGroups(std::vector<Group>& groupsToUpdate,
                const Protein<SasPdbAtom>& protein) const;

            private:
            std::vector<SerializableGroup> groups;
          };

          void makeGroups(float radius);
          void fillGroups(const std::string& sessionFileName, unsigned int timeStep);
          std::vector<Group> makeGroupsByDistance(const std::vector<Atom>& centers,
              float radius);
          std::vector<Group> makeGroupsByDistance(const std::vector<Atom>& centers,
              float radius,
              const std::vector<SasPdbAtom>& reference);
          Group makeGroupByDistance(const std::vector<Atom>& centers,
              const SasPdbAtom& atom, float radius);
          void pittpiRun();
          void clone(const Pittpi& pittpi); // Deprecated. Waiting for delegators...
#ifdef HAVE_PYMOD_SADIC
#if HAVE_PYMOD_SADIC == 1
          Protein<SasPdbAtom> runSadic(const Protein<SasPdbAtom>& structure) const;
#endif
#endif

          template<typename, typename> friend class Serializer;
          Gromacs gromacs;
          std::string sessionFileName;
          float radius;
          unsigned long threshold;
          Protein<SasPdbAtom> averageStructure;
          std::thread pittpiThread;
          mutable std::mutex statusMutex;
          mutable std::mutex nextStatusMutex;
          mutable std::condition_variable nextStatusCondition;
          mutable float __status;
          mutable std::string __statusDescription;
          bool sync;
          mutable std::mutex syncLock;
          bool abortFlag;
          std::vector<Pocket> pockets;
          std::vector<Group> groups;
        };

  template<typename Stream>
    void
    Pittpi::save(Stream& stream) const
    {
      Serializer<Stream> serializer(stream);
      serializer << radius;
      serializer << threshold;
      SerializableGroups serializableGroups(groups, averageStructure);
      serializer << serializableGroups;
      SerializablePockets serializablePockets(pockets, groups);
      serializer << serializablePockets;

      stream.close();
    }

  template<typename Stream>
    void
    Pittpi::load(Stream& stream)
    {
      assert(averageStructure.atoms().size() > 0);
      Serializer<Stream> serializer(stream);
      SerializableGroups serializableGroups;
      SerializablePockets serializablePockets;
      serializer >> radius;
      serializer >> threshold;
      serializer >> serializableGroups;
      serializer >> serializablePockets;

      serializableGroups.updateGroups(groups, averageStructure);
      serializablePockets.updatePockets(pockets, groups);

      sync = false;
    }

  template<class Serializer>
    void
    Pittpi::SerializablePockets::SerializablePocket::serialize(
        Serializer serializer)
    {
      serializer & groupIndex;
      serializer & startFrame;
      serializer & startPs;
      serializer & endFrame;
      serializer & endPs;
      serializer & width;
      serializer & openingFraction;
      serializer & averageNearFrame;
      serializer & averageNearPs;
      serializer & maxAreaFrame;
      serializer & maxAreaPs;
    }

  template<class Serializer>
    void
    Pittpi::SerializableGroups::SerializableGroup::serialize(
        Serializer serializer)
    {
      serializer & sas;
      serializer & zeros;
      serializer & referenceAtomIndex;
      serializer & referenceResIndex;
      serializer & residuesIndex;
    }

  template<typename Serializer>
    void
    Pittpi::SerializablePockets::serialize(Serializer serializer)
    {
      serializer & pockets;
    }

  template<typename Serializer>
    void
    Pittpi::SerializableGroups::serialize(Serializer serializer)
    {
      serializer & groups;
    }
}

#endif
