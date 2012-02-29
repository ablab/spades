/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2008) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////
/*
 * \file CreateLookupTab.cc
 * \author tsharpe
 * \date Jun 12, 2009
 *
 * \brief Writes a new-style lookup table.
 * A new-style lookup table is a feudal file of Locations (i.e., {contig#, offset} pairs).
 * The outer-vector index is the kmer, and the inner-vector contains the set of Locations for that kmer.
 *
 * This program just builds a singly-linked list of Locations for each kmer, and then walks each
 * list to write them out to a feudal file.
 */
#include "MainTools.h"
#include "lookup/LookupTabBuilder.h"

unsigned long gBaseCount;
unsigned int lineLen;

void progressReport( unsigned int contigNo, bvec const& bvec )
{
    unsigned long end = (gBaseCount + bvec.size())/1000000UL;
    for ( unsigned long iii = gBaseCount/1000000UL; iii < end; ++iii )
    {
        std::cout << '.';
        if ( ++lineLen == 60 )
        {
            std::cout << '\n';
            lineLen = 0;
        }
    }
    std::cout.flush();
    gBaseCount += bvec.size();
}

void processFile( char const* fileName, LookupTabBuilder& locList, bool is_fastb = false )
{
    std::cout << "Processing file " << fileName << std::endl;
    std::cout << "One dot per megabase processed:" << std::endl;

    size_t len = strlen(fileName);
    if ( is_fastb || len > 6 && !strcmp(fileName+len-6,".fastb") )
        locList.addFastb(fileName,progressReport);
    else
        locList.addFasta(fileName,progressReport);

    std::cout << std::endl;
}

int main( int argc, char** argv )
{
    RunTime();
    BeginCommandArguments;
    CommandArgument_UnsignedInt_OrDefault(K, 12);
    CommandArgument_String(SOURCE);
    CommandArgument_String(OUTPUT);

    // SOURCE is a fastb file, even if it does not end with ".fastb".
    CommandArgument_Bool_OrDefault( IS_FASTB, False );

    EndCommandArguments;

    LookupTabBuilder locList(K);
    if ( !SOURCE.EndsWith(".fof") )
        processFile( SOURCE.c_str(), locList, IS_FASTB );
    else // source file is a file of filenames, one per line
    {
        ifstream fof(SOURCE.c_str());
        char fileName[4096];
        while ( fof.getline(fileName,sizeof(fileName)) )
	  processFile( fileName, locList, IS_FASTB );
        fof.close();
    }

    locList.write(OUTPUT.c_str());
}
