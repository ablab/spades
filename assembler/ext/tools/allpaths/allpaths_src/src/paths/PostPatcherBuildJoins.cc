///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// MakeDepend: library OMP
// MakeDepend: cflags OMP_FLAGS

#include <omp.h>

#include "Basevector.h"
#include "MainTools.h"
#include "Qualvector.h"
#include "Superb.h"
#include "paths/UnipathFixerTools.h"
#include "feudal/BinaryStream.h"

int main( int argc, char *argv[] )
{
     RunTime( );

     BeginCommandArguments;
     CommandArgument_String(PRE);
     CommandArgument_String(DATA);
     CommandArgument_String(RUN);
     CommandArgument_String(SUBDIR);
     CommandArgument_Int(K);
     CommandArgument_String_OrDefault(SCAFFOLDS_IN, "linear_scaffolds");
     CommandArgument_String_OrDefault(JUMP_READS, "jump_reads_filt");
     CommandArgument_String_OrDefault(FRAG_READS, "frag_reads_filt");
     CommandArgument_String_OrDefault(POST_PATCH_DIR, "patch");
     CommandArgument_Bool_OrDefault(CHECKPOINT, False);
     EndCommandArguments;

     // Define directories and files.

     cout << Date( ) << ": entering PostPatcherBuildJoins" << endl;
     String data_dir = PRE + "/" + DATA, run_dir = PRE + "/" + DATA + "/" + RUN;
     String sub_dir = run_dir + "/ASSEMBLIES/" + SUBDIR;
     String ch_head 
          = sub_dir + "/" + POST_PATCH_DIR + "/PostPatcher." + SCAFFOLDS_IN + ".";
     String supers_file = sub_dir + "/" + SCAFFOLDS_IN + ".superb";
     String INGAP_file = ch_head + "INGAP", JINGAP_file = ch_head + "JINGAP";
     String reads_file = run_dir + "/" + FRAG_READS + ".fastb";
     String quals_file = run_dir + "/" + FRAG_READS + ".qualb";
     String jreads_file = run_dir + "/" + JUMP_READS + ".fastb";
     String jquals_file = run_dir + "/" + JUMP_READS + ".qualb";

     // Load data.

     cout << Date( ) << ": loading data" << endl;
     vec<superb> scaffolds;
     ReadSuperbs( supers_file, scaffolds );
     vec< vec< vec<opair> > > ingap, jingap;
     BinaryReader::readFile(INGAP_file.c_str(),&ingap);
     BinaryReader::readFile(JINGAP_file.c_str(),&jingap);
     String TIGS_file = ch_head + "TIGS";
     vecbasevector tigs(TIGS_file);

     // Build join data.

     cout << Date( ) << ": building join data" << endl;
     vec<int> gap_to_scaffold, gap_to_scaffold_pos;
     vec< vec<int> > gap_id( scaffolds.size( ) );
     uint64_t ngaps = 0;
     for ( int i = 0; i < scaffolds.isize( ); i++ )
     {    for ( int j = 0; j < scaffolds[i].Ngaps( ); j++ )
          {    gap_to_scaffold.push_back(i);
               gap_to_scaffold_pos.push_back(j);    
               gap_id[i].push_back(ngaps++);    }    }

     cout << Date() << ": looking for joins across " << ngaps << " scaffold gaps" << endl;

     size_t index_counter = 0;
     vec< vec<size_t> > index(ngaps);
     vecbasevector extenders_b;
     vecqualvector extenders_q;
     vec< vec<opair> > extenders(ngaps);

     cout << Date() << ": processing fragment reads..." << endl;
     {
       cout << Date( ) << ":   loading reads.         mem =  " << ToStringAddCommas( MemUsageBytes( ) ) << endl;
       vecbasevector reads(reads_file);
       cout << Date( ) << ":   building extenders.    mem =  " << ToStringAddCommas( MemUsageBytes( ) ) << endl;
       for ( int i = 0; i < scaffolds.isize( ); i++ ) {
	 for ( int j = 0; j < scaffolds[i].Ngaps( ); j++ ) {
	   UniqueSort( ingap[i][j] );
	   for ( int k = 0; k < ingap[i][j].isize( ); k++ ) {
	     const opair& p = ingap[i][j][k];
	     basevector b = reads[p.id2];
	     b.ReverseComplement( );
	     int g = gap_id[i][j];
	     index[g].push_back(index_counter, index_counter + 1);
	     index_counter+=2;
	     extenders_b.push_back( reads[p.id1] );
	     extenders_b.push_back( b ); 
	   } 
	 }
       }
       cout << Date( ) << ":   found " << extenders_b.size() << " fragment extenders" << endl; 
       cout << Date( ) << ":   releasing reads.       mem =  " << ToStringAddCommas( MemUsageBytes( ) ) << endl;
     }

     cout << Date() << ": processing fragment quals..." << endl;
     {
       cout << Date( ) << ":   loading quals.         mem =  " << ToStringAddCommas( MemUsageBytes( ) ) << endl;
       vecqualvector quals(quals_file);
       cout << Date( ) << ":   building extenders.    mem =  " << ToStringAddCommas( MemUsageBytes( ) ) << endl;
       for ( int i = 0; i < scaffolds.isize( ); i++ ) {
	 for ( int j = 0; j < scaffolds[i].Ngaps( ); j++ ) {
	   for ( int k = 0; k < ingap[i][j].isize( ); k++ ) {
	     const opair& p = ingap[i][j][k];
	     int g = gap_id[i][j];
	     extenders_q.push_back( quals[p.id1] );
	     extenders_q.push_back( Reverse( quals[p.id2] ) );    
	     extenders[g].push_back(p);
	   }
	 }
       }
       cout << Date( ) << ":   releasing quals.       mem =  " << ToStringAddCommas( MemUsageBytes( ) ) << endl;
     }

     Destroy(ingap);

     size_t frag_extender_count = extenders_b.size();
     cout << Date() << ": processing jumping reads..." << endl;
     {
       cout << Date( ) << ":   loading reads.         mem =  " << ToStringAddCommas( MemUsageBytes( ) ) << endl;
       vecbasevector jreads(jreads_file);
       cout << Date( ) << ":   building extenders.    mem =  " << ToStringAddCommas( MemUsageBytes( ) ) << endl;
       for ( size_t i = 0; i < jreads.size( ); i++ )
	 jreads[i].ReverseComplement( );
       for ( int i = 0; i < scaffolds.isize( ); i++ ) {
	 for ( int j = 0; j < scaffolds[i].Ngaps( ); j++ ) {
           UniqueSort( jingap[i][j] ); 
	   for ( int k = 0; k < jingap[i][j].isize( ); k++ ) {
	     const opair& p = jingap[i][j][k];
	     basevector b = jreads[p.id2];
	     b.ReverseComplement( );
	     int g = gap_id[i][j];
	     index[g].push_back(index_counter, index_counter + 1);
	     index_counter+=2;
	     extenders_b.push_back( jreads[p.id1] );
	     extenders_b.push_back( b );
	   }
	 }
       }
       cout << Date( ) << ":   found " << extenders_b.size() - frag_extender_count << " jumping extenders" << endl; 
       cout << Date( ) << ":   releasing reads.       mem =  " << ToStringAddCommas( MemUsageBytes( ) ) << endl;
     }

     cout << Date() << ": processing jumping quals..." << endl;
     {
       cout << Date( ) << ":   loading quals.         mem =  " << ToStringAddCommas( MemUsageBytes( ) ) << endl;
       vecqualvector jquals(jquals_file);
       cout << Date( ) << ":   building extenders.    mem =  " << ToStringAddCommas( MemUsageBytes( ) ) << endl;
       for ( size_t i = 0; i < jquals.size( ); i++ )
	 jquals[i].ReverseMe( );
       cout << Date( ) << ":   finished reversing.    mem =  " << ToStringAddCommas( MemUsageBytes( ) ) << endl;
       for ( int i = 0; i < scaffolds.isize( ); i++ ) {
	 for ( int j = 0; j < scaffolds[i].Ngaps( ); j++ ) {
	   for ( int k = 0; k < jingap[i][j].isize( ); k++ ) {
	     const opair& p = jingap[i][j][k];
	     int g = gap_id[i][j];
	     extenders_q.push_back( jquals[p.id1] );
	     extenders_q.push_back( Reverse( jquals[p.id2] ) );    
	     extenders[g].push_back(p);
	   }
	 }
       }
       cout << Date( ) << ":   releasing quals.       mem =  " << ToStringAddCommas( MemUsageBytes( ) ) << endl;
     }
     cout << Date( ) << ": all extensions built.   mem =  " << ToStringAddCommas( MemUsageBytes( ) ) << endl;

     Destroy(jingap);

     ForceAssert(index_counter == extenders_b.size());

     cout << Date( ) << ": writing temp extenders" << endl;

     if (CHECKPOINT) {
       extenders_b.WriteAll(ch_head + "e_b");
       extenders_q.WriteAll(ch_head + "q_b");
       String index_file = ch_head + "index";
       BinaryWriter::writeFile( index_file.c_str(), index );
       String e_p_file = ch_head + "e_p";
       BinaryWriter::writeFile( e_p_file.c_str(), extenders );
     }

     // Package input data for joins.

     cout << Date( ) << ": package input data for joins" << endl;

     vec<size_t> join_data_offsets;
     join_data_offsets.reserve(ngaps);

     String JOINDATA_file = ch_head + "JOINDATA", JDLEN_file = ch_head + "JDLEN";
     BinaryWriter jdWriter(JOINDATA_file.c_str(),false);
     jdWriter.write(PCottageJoinData::HEADER);
     size_t offset = sizeof(PCottageJoinData::HEADER);

     PCottageJoinData joinData;
     for ( uint64_t i = 0; i < ngaps; i++ )
     {
       const superb& S = scaffolds[ gap_to_scaffold[i] ];
       int p = gap_to_scaffold_pos[i];
       joinData.sep = S.Gap(p);
       joinData.dev = S.Dev(p);
       joinData.L = tigs[S.Tig(p)];
       joinData.R = tigs[S.Tig(p+1)];

       // Combine the fragment and jump data.

       vec<size_t> const& idxEntry = index[i];
       joinData.reads.clear().reserve(idxEntry.size());
       joinData.quals.clear().reserve(idxEntry.size());
       typedef vec<size_t>::const_iterator Itr;
       for ( Itr itr(idxEntry.begin()), end(idxEntry.end()); itr != end; ++itr )
       {
	   joinData.reads.push_back(extenders_b[*itr]);
	   joinData.quals.push_back(extenders_q[*itr]);
       }

       vec<opair> const& extEntry = extenders[i];
       joinData.pairs.clear();
       joinData.pairs.reserve(extEntry.size());
       typedef vec<opair>::const_iterator XItr;
       for ( XItr itr(extEntry.begin()), end(extEntry.end()); itr!=end; ++itr )
       {
           opair const& p = *itr;
           joinData.pairs.push_back(make_pair<int,int>(p.sep,p.dev));
       }

       join_data_offsets.push_back(offset);
       offset += jdWriter.write(joinData);
     }
     jdWriter.close();
     BinaryWriter::writeFile(JDLEN_file.c_str(),join_data_offsets);

     cout << Date( ) << ": join data completed" << endl;

}
