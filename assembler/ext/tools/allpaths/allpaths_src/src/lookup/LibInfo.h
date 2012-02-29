///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////
/*
 * \file LibInfo.h
 * \author tsharpe
 * \date Mar 16, 2010
 *
 * \brief
 */
#ifndef LIBINFO_H_
#define LIBINFO_H_

#include "String.h"
#include "TokenizeString.h"
#include "Vec.h"
#include "system/System.h"
#include <map>
#include <fstream>

struct LibInfo
{
    LibInfo() : mMean(0), mStdDev(0) {}
    LibInfo( int mean, int stdDev ) : mMean(mean), mStdDev(stdDev) {}

    // compiler-suppied copying and destructor are OK

    int mMean;
    int mStdDev;
};

class LibInfoDB
{
public:
    LibInfoDB( String const& libInfoFilename )
    {
        ifstream in(libInfoFilename.c_str());
        char buf[8192];
        vec<String> toks;
        toks.reserve(3);
        while ( in.getline(buf,sizeof(buf)) )
        {
            unsigned len = strlen(buf);
            if ( len && buf[0] != '#' )
            {
                int nToks = Tokenize(buf,toks);
		// Trim off anything after a '#' token.
		for ( int i = 0; i < nToks; i++ )
		  if ( toks[i] == "#" ) {
		    toks.resize(i);
		    nToks = i;
		    break;
		  }
		
                long mean;
                long stdDev;
                if ( nToks && nToks != 3 ||
                        !toks[1].IsInt(&mean) ||
                        !toks[2].IsInt(&stdDev) )
                    FatalErr("Couldn't interpret this gibberish from the libInfo file, "
                             << libInfoFilename << ": " << buf );
                mMap[toks[0]] = LibInfo(mean,stdDev);
            }
        }
        if ( !in.eof() )
            FatalErr("Couldn't read libInfo file: " << libInfoFilename );
        in.close();
    }

    // compiler-supplied copying and destructor are OK

    LibInfo const* getInfo( String const& libName ) const
    { std::map<String,LibInfo>::const_iterator itr = mMap.find(libName);
      LibInfo const* result = 0;
      if ( itr != mMap.end() ) result = &itr->second;
      return result; }

    String getLibraryName(int mean, int stdDev) const {
      for (std::map<String,LibInfo>::const_iterator itr = mMap.begin(); 
	   itr != mMap.end(); 
	   ++itr) {
	if ((itr->second).mMean == mean && (itr->second).mStdDev == stdDev) return itr->first;
      }
      return String();
    }
  
private:
    std::map<String,LibInfo> mMap;
};

#endif /* LIBINFO_H_ */
