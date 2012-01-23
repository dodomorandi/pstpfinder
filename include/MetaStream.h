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
#include <iosfwd>
#include <boost/iostreams/categories.hpp>

using namespace std;
namespace io = boost::iostreams;

namespace PstpFinder
{

  class MetaStream
  {
    public:
      typedef char char_type;
      typedef io::source_tag category;

      MetaStream();
      MetaStream(ifstream& modifiableStream, streampos begin = -1,
                 streampos end = -1);
      MetaStream(const MetaStream& metaStream);

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

      streamsize read(char* data, streamsize length);
      MetaStream& seekg(streamsize pos);
      MetaStream& seekg(streamsize off, ios_base::seek_dir way);

      streamsize tellg() const;
      bool eof() const;

    private:
      ifstream& inputStream;
      ifstream nullStream;
      streampos streamBegin;
      streampos streamEnd;
      const bool valid;
      bool eofTrigger;

      template<typename T> void getFromStream(T& out);
  };

} /* namespace PstpFinder */
#endif /* METASTREAM_H_ */
