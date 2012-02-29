///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////
/*
 * \file PathReads.cc
 * \author tsharpe
 * \date Sep 13, 2010
 *
 * \brief
 */
#include "kmers/ReadPatherDefs.h"
#include "MainTools.h"

int main( int argc, char** argv )
{
    RunTime();
    BeginCommandArguments;
    CommandArgument_String_OrDefault(PRE,".");
    CommandArgument_String_OrDefault(DATA,".");
    CommandArgument_String_OrDefault(RUN,".");
    CommandArgument_String_Doc(READS_IN,"Base name of fastb file.");
    CommandArgument_String_OrDefault_Doc(GRAPHINFO,"","KmerDB to path on.");
    CommandArgument_UnsignedInt_OrDefault_Doc(COVERAGE_ESTIMATE,50,
            "Used to size hash table.  Assume 50x coverage, in the absence of "
            "a better guess.");
    CommandArgument_UnsignedInt_OrDefault(NUM_THREADS,0);
    CommandArgument_Bool_OrDefault(VALIDATE,True);
    CommandArgument_Bool_OrDefault(WRITE_KMER_COUNTS,False);
    CommandArgument_LongLong_OrDefault(MAX_MEMORY_GB,0);
    EndCommandArguments;

    unsigned const K = 28;

    SetMaxMemory(MAX_MEMORY_GB<<30);

    String readsFile = PRE + '/' + DATA + '/' + RUN + '/' + READS_IN + ".fastb";
    if ( GRAPHINFO == "" )
        GRAPHINFO = UnipathGraph<K>::getInfoFilename(readsFile);

    if ( NUM_THREADS < 1 || NUM_THREADS > processorsOnline() )
        NUM_THREADS = processorsOnline();

    if ( IsRegularFile(GRAPHINFO) )
        PathCollection<K>::create(readsFile,GRAPHINFO,VALIDATE,NUM_THREADS);
    else
    {
        size_t nKmers = 4 * MasterVec<bvec>::MastervecFileRawCount(readsFile) /
                            COVERAGE_ESTIMATE;
        PathCollection<K>::create(readsFile,VALIDATE,NUM_THREADS,nKmers,
                                WRITE_KMER_COUNTS);
    }

    PathCollection<K> pc(PathCollection<K>::getInfoFilename(readsFile),GRAPHINFO);
    if ( VALIDATE )
        pc.validate(readsFile);

    std::cout << Date() << ": Done." << std::endl;
}
