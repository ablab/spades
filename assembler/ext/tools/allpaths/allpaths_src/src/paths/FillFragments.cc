///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////
/*
 * \file FillFragments.cc
 * \author tsharpe
 * \date Sep 14, 2010
 *
 * \brief
 */
#include "MainTools.h"
#include "feudal/BitVec.h"
#include "paths/FragmentFillerDefs.h"
#include "util/SearchFastb2Core.h"
#include <algorithm>
#include <fstream>

#include "reporting/PerfStat.h"

int main( int argc, char** argv )
{
    String empty;
    RunTime();
    BeginCommandArguments;
    CommandArgument_String_OrDefault(PRE,".");
    CommandArgument_String_OrDefault(DATA,".");
    CommandArgument_String_OrDefault(RUN,".");
    CommandArgument_UnsignedInt_OrDefault_Doc(NUM_THREADS, 0, 
            "Number of threads to use (use all available processors if set to 0)");
    CommandArgument_String_OrDefault(READS_IN,"frag_reads_corr");
    CommandArgument_String_OrDefault(PAIRS_IN,empty);
    CommandArgument_String_OrDefault(READS_OUT,"filled_reads");
    CommandArgument_String_OrDefault(PAIRS_OUT,empty);
    CommandArgument_String_OrDefault(PAIR_IDS,empty);
    CommandArgument_UnsignedInt_OrDefault_Doc(MIN_OVERLAP,0,
            "No. of bases overlap in valid glue.");
    CommandArgument_UnsignedInt_OrDefault_Doc( MAX_CLIQ, 64000,
            "Maximum amount of evidence to wade through.");
    CommandArgument_Double_OrDefault_Doc(MAX_STRETCH,3.,
            "Maximum number of SDs in a valid pair separation.");
    CommandArgument_Bool_OrDefault_Doc(UNIQUE_ONLY,False,
            "Emit closures only when unique.");
    CommandArgument_Bool_OrDefault_Doc(PRECORRECT_LIBSTATS,False,
            "Pre-sample input to estimate library stats.");
    CommandArgument_String_OrDefault(GRAPHINFO,empty);
    CommandArgument_String_OrDefault(PATHINFO,empty);
    CommandArgument_UnsignedInt_OrDefault_Doc( MAX_GLUE, 100,
            "Maximum number of glue reads to consider.");
    CommandArgument_String_OrDefault_Doc(EXTRA_FILLS, empty,
            "extra fastb file of bogus filled fragments to be appended to the "
            "file of real filled fragments - temporary hack");
    CommandArgument_String_OrDefault_Doc(XREF_FASTB, empty,
            "Fastb file containing extended reference genome.  If defaulted "
            "no accuracy check is done.");
    EndCommandArguments;

    unsigned const K = 28;

    // Thread control
   
    NUM_THREADS = configNumThreads(NUM_THREADS);
 
    String runDir = PRE + '/' + DATA + '/' + RUN + '/';
    String readsFile = runDir + READS_IN + ".fastb";
    String filledFile = runDir + READS_OUT + ".fastb";

    if ( MIN_OVERLAP < K )
        MIN_OVERLAP = K;
    long kOverlap = MIN_OVERLAP - K + 1;

    if ( PAIRS_IN == empty )
        PAIRS_IN = READS_IN;
    PAIRS_IN = runDir + PAIRS_IN + ".pairs";
    if ( PAIRS_OUT == empty )
        PAIRS_OUT = READS_IN + "_cpd";
    PAIRS_OUT = runDir + PAIRS_OUT + ".pairs";

    if ( GRAPHINFO == empty )
        GRAPHINFO = UnipathGraph<K>::getInfoFilename(readsFile);
    if ( PATHINFO == empty )
        PATHINFO = PathCollection<K>::getInfoFilename(readsFile);

    std::cout << Date() << " Loading paths." << std::endl;
    FragmentFiller<K> ff(PATHINFO, GRAPHINFO, kOverlap, MAX_CLIQ, MAX_GLUE,
            !UNIQUE_ONLY);
    std::cout << Date() << " Loading pairs." << std::endl;
    PairsManager pm(PAIRS_IN);
    if ( PRECORRECT_LIBSTATS )
    {
        LibSepsColl lsc(pm, K, MAX_STRETCH);
        size_t nPairs = std::min(pm.nPairs(), 10000ul);
        for ( size_t idx = 0; idx < nPairs; ++idx )
            ff.processPair(lsc, idx, 0);
        if ( lsc.getNFills() * 10 < nPairs )
            std::cout
                    << "No library parameter pre-correction: too few pairs closed."
                    << std::endl;
        else
        {
            lsc.updatePM();
            std::cout << "Results from library parameter pre-correction:"
                    << std::endl;
            pm.printLibraryStats(std::cout);
        }
    }
    LibSepsColl lsc(pm, K, MAX_STRETCH);
    size_t nFragsFilled = 0;

    if ( PAIR_IDS != empty )
    {
        vec<int> ids;
        ParseIntSet(PAIR_IDS, ids);

        FragmentFiller<K>::OutVec out;
        IncrementalWriter<bvec> writer(filledFile.c_str(), ids.size());
        ofstream infoOS(filledFile.ReplaceExtension(".fastb", ".info").c_str());
        typedef vec<int>::iterator Itr;
        for ( Itr itr(ids.begin()), end(ids.end()); itr != end; ++itr )
        {
            size_t pairId = *itr;
            if ( pairId >= pm.nPairs() )
                std::cout << "Ignoring illegal pairId " << pairId << std::endl;
            else
            {
                ff.processPair(lsc, *itr, &out);
                typedef FragmentFiller<K>::OutVec::iterator OItr;
                OItr oEnd(out.end());
                for ( OItr oItr(out.begin()); oItr != oEnd; ++oItr )
                {
                    writer.add(oItr->getFill());
                    oItr->write(infoOS,nFragsFilled++) << '\n';
                    out.clear();
                }
            }
        }
        writer.close();
        infoOS.close();
    }
    else
    {
        nFragsFilled = ff.processPairs(lsc,filledFile,NUM_THREADS);
        if ( nFragsFilled * 10 < pm.nPairs() )
	{
	  std::cout << "No library parameter adjustment:  too few pairs closed."
		    << std::endl;
	  FatalErr("Less than 10% of fragment pairs were filled.\n"
		   "There may be a problem with the library.");
	}
        else
        {
            lsc.updatePM();
            pm.printLibraryStats(std::cout);
        }
        pm.Write(PAIRS_OUT);
    }

    PerfStat::log() << std::fixed << std::setprecision(1)
                    << PerfStat("frac_filled_pairs",
                                "% of fragment pairs that were filled",
                                100.*float(nFragsFilled)/float(pm.nPairs()));


    if ( XREF_FASTB != empty )
    {
        BitVec isMatch;
        SearchFastb2(filledFile,XREF_FASTB,96,0,&isMatch,0,.90,False);
        PerfStat::log() << std::fixed << std::setprecision(1)
                        << PerfStat("frac_accurate_pairs",
                                    "% of fills that match genome",
                                    100.*isMatch.Sum()/isMatch.size());
    }

    if ( EXTRA_FILLS != empty )
    {   vecbasevector fills(filledFile), extras(EXTRA_FILLS);
        fills.Append(extras);
        fills.WriteAll(filledFile);
        nFragsFilled += extras.size( );}

    PairsManager pmEmpty(nFragsFilled);
    pmEmpty.Write(filledFile.ReplaceExtension(".fastb",".pairs"));

    std::cout << Date() << " Done." << std::endl;
}
