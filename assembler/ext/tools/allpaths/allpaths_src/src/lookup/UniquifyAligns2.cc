// Copyright (c) 2000-2003 Whitehead Institute for Biomedical Research

/**
   Program: UniquifyAligns2
   
   Given alignments of reads to assembly, try to keep only the correct alignments
   and throw out incorrect ones; also, identify those alignments where we're highly
   sure that the alignment is unique and correct.

   Primarily, we use pairing information.  

   We look for situations where there is a read pair, neither read of which is 
   placed uniquely, but such that there is only one placement of each read 
   consistent with the pairing.
   
   Produce as output a list of all the unique alignments, ordered by read id, 
   including those disambiguated in the ways just described.

   To accept a read pair alignment as unique, it must be within 
   MAX_DISCREP_TO_ACCEPT standard deviations of the predicted value, and there must
   not be any other pair alignments within MAX_DISCREP_TO_ACCEPT standard 
   deviations of the predicted value.

   See also: UniquifyAligns

   Inputs:

   reads.fastb
   reads.pairs
*/


#include <strstream>

#include "MainTools.h"
#include "FastIfstream.h"
#include "math/Functions.h"
#include "PackAlign.h"
#include "PairsManager.h"
#include "CommonSemanticTypes.h"

#include "lookup/LookAlign.h"
#include "lookup/LookupTable.h"


class placement {

     public:

     int ind1, ind2;
     float discrep;

     placement( ) { }
     placement( int ind1_arg,
		int ind2_arg,
		float discrep_arg )
       : ind1(ind1_arg),
	 ind2(ind2_arg),
	 discrep(discrep_arg)
     { }

};

struct placement_cmp_by_discrep {
public:
  bool operator() ( const placement& p1, const placement& p2 ) const {
    return p1.discrep < p2.discrep;
  }
};


int main( int argc, char *argv[] )
{
     RunTime( );

     BeginCommandArguments;
     CommandArgument_String(PRE);
     CommandArgument_String(DATA);
     CommandArgument_String_OrDefault(RUN, "");
     CommandArgument_Bool_OrDefault(DATA_IN_RUN_DIR, False);
     CommandArgument_Int_OrDefault_Doc(MAX_DISCREP_TO_ACCEPT, 3,
				       "max # of standard deviations that a read separation "
				       "may have from the mean, for us to accept a pair placement "
				       "on the reference as plausible." );
     CommandArgument_Int_OrDefault_Doc(MAX_DISCREP_TO_CONSIDER, 4,
				       "there must not be any other pair placements with "
				       "separations within this # of stddevs of the mean "
				       "for us to accept a pair placement.");
     CommandArgument_String_OrDefault_Doc( READS, "reads",
					   "prefix of the file containing the reads that were "
					   "aligned to the reference.  normally the file is reads.fastb "
					   "in the run directory." );
     CommandArgument_String(ALL_ALIGNS_IN);
     CommandArgument_String(UNIQ_ALIGNS_OUT);
     CommandArgument_Bool_OrDefault(VERBOSE, False);
     EndCommandArguments;

     // Load some assembly data.

     String data_dir = PRE + "/" + DATA;
     String run_dir = data_dir + "/" + RUN;

     String working_dir = (DATA_IN_RUN_DIR ? run_dir : data_dir);

     cout << Date() << " Loading pair information..." << endl;
     int nreads = MastervecFileObjectCount( working_dir + "/" + READS + ".fastb" );
     
     PairsManager pairs( working_dir + "/" + READS + ".pairs" );
     cout << Date() << " Loaded " << pairs.nPairs() << " pairs." << endl;

     // Possibly here: remove aligns where the combination
     // of errors implied by the align is highly unlikely
     // (based on quality scores / trustedness / intensisites / etc).

     // Build an index: for each read, aligns of that read.

     vec<look_align> aligns;
     vec< vec<align_id_t> > aligns_index;
     filename_t alignsFile = working_dir + "/" + ALL_ALIGNS_IN;

     cout << Date() << " Loading aligns for " << nreads << " reads..." << endl;
     LoadLookAligns( alignsFile, aligns, aligns_index, nreads );
     cout << Date() << " Loaded " << aligns.isize() << " aligns." << endl;

     vec< look_align > uniqueAligns;

     // Keep all aligns of unpaired reads
     for ( longlong readId = 0 ; readId < nreads; readId++ ) 
       if ( pairs.isUnpaired( readId ) &&
	    aligns_index[ readId ].solo() ) 
	 uniqueAligns.push_back( aligns[ aligns_index[ readId ].front() ] );

     // Disambiguate pairs.

     cout << " Disambiguating " << pairs.nPairs() << " pairs..." << endl;
     int numDisambiguated = 0;
     for ( size_t pairId = 0; pairId < pairs.nPairs( ); pairId++ ) {
       if ( !(pairId % 100000) ) {
	 cout << Date() << " i=" << pairId << " disamb=" << numDisambiguated << endl;
       }

       longlong id1 = pairs.ID1( pairId ), id2 = pairs.ID2( pairId );
       int sep = pairs.sep( pairId ), sd = pairs.sd( pairId );
       
       int min_sep_accept = max( 0, sep - MAX_DISCREP_TO_ACCEPT * sd ),
	 max_sep_accept = sep + MAX_DISCREP_TO_ACCEPT * sd;

       int min_sep_consider = max( 0, sep - MAX_DISCREP_TO_CONSIDER * sd ),
	 max_sep_consider = sep + MAX_DISCREP_TO_CONSIDER * sd;

       triple< align_id_t, align_id_t, nbases_t > consideredPlacement( -1, -1, 0 );
       int numConsideredPlacements = 0;
       
       for ( int j1 = 0; j1 < aligns_index[id1].isize( ) && numConsideredPlacements <= 1; j1++ ) {
	 align_id_t ind1 = aligns_index[id1][j1];
	 const look_align& la1 = aligns[ind1];
	 for ( int j2 = 0; j2 < aligns_index[id2].isize( ) && numConsideredPlacements <= 1; j2++ ) {
	   align_id_t ind2 = aligns_index[id2][j2];
	   const look_align& la2 = aligns[ind2];
	   
	   if ( la1.target_id != la2.target_id ) continue;
	   if ( la1.rc1 == la2.rc1 ) continue;
	       
	   nbases_t sep = 
	     !la1.rc1 ? la2.a.pos2( ) - la1.a.Pos2( ) : la1.a.pos2( ) - la2.a.Pos2( );
	   
	   if ( min_sep_consider < sep  && sep < max_sep_consider ) {
	     consideredPlacement = make_triple( ind1, ind2, sep );
	     numConsideredPlacements++;
	   }
	 }
       }  
       if ( numConsideredPlacements == 1 &&
	    min_sep_accept < consideredPlacement.third && consideredPlacement.third < max_sep_accept ) {
	      uniqueAligns.push_back( aligns[ consideredPlacement.first ] );
	      uniqueAligns.push_back( aligns[ consideredPlacement.second ] );
       }
     }

     // Output unique alignments.

     cout << Date() << " Writing " << uniqueAligns.size() << " unique aligns..." << endl;
     WriteLookAligns( working_dir + "/" + UNIQ_ALIGNS_OUT  , uniqueAligns );
}
