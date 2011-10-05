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

#ifndef METASTREAM_H_
#define METASTREAM_H_

#include <string>
#include <fstream>

using namespace std;

namespace Gromacs
{

  class MetaStream
  {
    public:
      MetaStream(const string& fileName, streampos begin = 0,
                 streampos end = -1);

      MetaStream& operator >>(string& out);
      MetaStream& operator >>(bool& out);
      MetaStream& operator >>(char& out);
      MetaStream& operator >>(unsigned char& out);
      MetaStream& operator >>(short& out);
      MetaStream& operator >>(unsigned short& out);
      MetaStream& operator >>(int& out);
      MetaStream& operator >>(unsigned int& out);
      MetaStream& operator >>(long& out);
      MetaStream& operator >>(unsigned long& out);
      MetaStream& operator >>(void*& out);
      MetaStream& operator >>(float& out);
      MetaStream& operator >>(double& out);
      MetaStream& operator >>(long double& out);

      MetaStream& seekg(long pos);
      MetaStream& seekg(long off, ios_base::seek_dir way);

      long tellg();
      bool eof() const;

    private:
      ifstream inputStream;
      streampos streamBegin;
      streampos streamEnd;
      streampos currentPosition;

      void checkEof() const;
      template<typename T> void getFromStream(T& out);
  };

} /* namespace Gromacs */
#endif /* METASTREAM_H_ */
