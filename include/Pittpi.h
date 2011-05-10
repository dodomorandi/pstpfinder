#ifndef _PITTPI_H
#define _PITTPI_H

#include "Gromacs.h"
#include <vector>

namespace Gromacs
{
  /**
   * @brief Related to the method developed by Matteo De Chiara
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
    Pittpi(const Gromacs& gromacs, float radius, unsigned long threshold);
  private:
    std::vector<const Residue*> makeGroups(unsigned long threshold);

    Protein averageStructure;
  };
}

#endif
