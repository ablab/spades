///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2011) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// MakeDepend: library OMP
// MakeDepend: cflags OMP_FLAGS
#include <omp.h>

#include "CoreTools.h"
#include "VecUtilities.h"
#include "ParallelVecUtilities.h"
#include "efasta/EfastaTools.h"
#include "paths/AssemblyCleanupTools.h"


Bool scompare( const superb& s1, const superb& s2 ){
  int l1 = s1.FullLength();
  int l2 = s2.FullLength();
  return (l1 > l2);
}


Assembly::Assembly( const String in_superb_file, const String in_efasta_file ){
  // Loading scaffolds
  cout << Date( ) << ": loading superb file" << endl;
  ReadSuperbs( in_superb_file, scaffolds );
  
  // reading contig information
  LoadEfastaIntoStrings(in_efasta_file, efastas);
  tigMap.resize( efastas.size() );
  for ( int tid = 0; tid < efastas.isize(); tid++ )
    tigMap[tid] = ToString(tid);

  scaffMap.resize( scaffolds.size() );
  for ( int sid = 0; sid < scaffolds.isize(); sid++ )
    scaffMap[sid] = ToString(sid);
}


Assembly::Assembly( const vec<superb>& scaffoldsIn, const vec<efasta>& efastasIn ){

  scaffolds = scaffoldsIn;
  efastas   = efastasIn;
  tigMap.resize( efastas.size() );
  for ( int tid = 0; tid < efastas.isize(); tid++ )
    tigMap[tid] = ToString(tid);

  scaffMap.resize( scaffolds.size() );
  for ( int sid = 0; sid < scaffolds.isize(); sid++ )
    scaffMap[sid] = ToString(sid);
}


size_t Assembly::scaffoldsTotLen() const{
  size_t len = 0;
  for ( size_t is = 0; is < scaffolds.size(); is++ )
    len += scaffolds[is].FullLength();
  return len;
}

size_t Assembly::scaffoldsRedLen() const{
  size_t len = 0;
  for ( size_t is = 0; is < scaffolds.size(); is++ )
    len += scaffolds[is].ReducedLength();
  return len;
}

size_t Assembly::scaffoldsNtigs() const{
  size_t ntigs = 0;
  for ( size_t is = 0; is < scaffolds.size(); is++ )
    ntigs += scaffolds[is].Ntigs();
  return ntigs;
}

// check integrity of scafolds and contigs data: contig size in superb == contig size in efasta,
//  each contig used once and only once
void Assembly::check_integrity() const{
  
  cout << Date() << ": checking integrity" << endl;
  vec<int> tigs_used( efastas.size(), 0);
  for ( size_t i = 0; i < efastas.size( ); i++ )
  {    vec<String> s(1);
       s[0] = efastas[i];
       ValidateEfastaRecord(s);    }
  for ( size_t si = 0; si < scaffolds.size(); si++ ){
    const superb & s = scaffolds[si];
    ForceAssertGt( s.Ntigs(), 0 );
    for ( int tpos = 0; tpos < s.Ntigs(); tpos++ ){
      size_t tid = s.Tig(tpos);
      ForceAssertLt( tid, efastas.size() );
      if ( efastas[tid].Length1() != s.Len(tpos) ){
	PRINT5( si, tpos, tid, s.Len(tpos), efastas[tid].Length1() );
	ForceAssertEq( efastas[tid].Length1(), s.Len(tpos) );
      }
      tigs_used[tid]++;
    }
  }
  vec<size_t> unused_tigs, overused_tigs;
  for ( size_t tid = 0; tid < tigs_used.size(); tid++ ){
    if ( tigs_used[tid] == 0 )
      unused_tigs.push_back( tid );
    else if ( tigs_used[tid] > 1 )
      overused_tigs.push_back(tid);
    
  }
  
  if ( unused_tigs.size() > 0 || overused_tigs.size() > 0 ){
    
    if ( unused_tigs.size() > 0 ){
      int max_un_size = efastas.at( unused_tigs[0] ).Length1();
      for ( size_t i = 0; i < unused_tigs.size(); i++ )
	if (  efastas.at( unused_tigs[i] ).Length1() > max_un_size )
	  max_un_size = efastas.at( unused_tigs[i] ).Length1();

      cout << "maximum size of unused contig is : " << max_un_size << endl;
    }

    PRINT2( unused_tigs.size(), overused_tigs.size() );
    ForceAssert( unused_tigs.size() == 0 && overused_tigs.size() == 0 );
  }
  
  return;
}


void Assembly::remove_small_scaffolds(const int min_scaffold_size) {
  cout << Date() << " removing small scaffolds: " << endl;
  cout << "initial number of scaffolds = " << scaffolds.size() << endl;
  for ( int si = 0; si < scaffolds.isize(); si++ )
    if ( scaffolds.at(si).ReducedLength() < min_scaffold_size ){
      scaffolds.erase( scaffolds.begin() + si );
      scaffMap.erase( scaffMap.begin() + si );
      si--;      
    }
  for ( int si = 0; si < scaffolds.isize(); si++ )
    ForceAssertGe( scaffolds.at(si).ReducedLength(), min_scaffold_size );
  cout << "final number of scaffolds = " << scaffolds.size() << endl;
}

void Assembly::remove_contigs( const vec<Bool>& to_remove )
{
  vec<int> offsets( efastas.size(), 0 );

  int offset = 0;
  for ( size_t tid = 0; tid < efastas.size(); tid++ ){
    if ( to_remove[tid] ){
      offsets[tid] = -1;
      offset++;
    }
    else{ offsets[tid] = offset; }
  }
    
  ForceAssertEq( offsets.size(), efastas.size() );
  for ( int tid = 0; tid < offsets.isize(); tid++ ){
    if ( offsets[tid] > 0 ){
      efastas[ tid - offsets[tid] ] = efastas[tid];
      tigMap[tid - offsets[tid] ] = tigMap[tid];
    }      
  }
  efastas.resize( efastas.size() - offset );
  tigMap.resize( tigMap.size() - offset );
  ForceAssertEq( efastas.size(), tigMap.size() );
  
  for ( size_t si = 0; si < scaffolds.size(); si++ ){
    for ( int tpos = 0; tpos < scaffolds[si].Ntigs(); tpos++ ){
      int tid = scaffolds[si].Tig(tpos);
      if ( offsets[tid] >= 0 ){
	int newtid = tid - offsets[tid];
	ForceAssertGe( newtid, 0 );
	scaffolds[si].SetTig( tpos, newtid );
      }
      else{
	scaffolds[si].RemoveTigByPos( tpos );
	tpos--;
      }
    }
  }
}

   
void Assembly::remove_small_contigs( const int min_size_solo, const int min_size_in ){
  cout << Date() << " removing small contigs and renumbering" << endl;
  vec<int> offsets( efastas.size(), 0 );
  vec<Bool> tigsToRemove( efastas.size(), False);
  for ( size_t si = 0; si < scaffolds.size(); si++ ){
    if ( scaffolds[si].Ntigs() == 1 ){
      if ( scaffolds[si].Len(0) < min_size_solo )
	tigsToRemove[ scaffolds[si].Tig(0) ] = True;
    }
    else{
      for ( int tpos = 0; tpos < scaffolds[si].Ntigs(); tpos++ ){
	if ( scaffolds[si].Len(tpos) < min_size_in )
	  tigsToRemove[ scaffolds[si].Tig(tpos) ] = True;
      }
    }
  }
  remove_contigs(tigsToRemove);
}

void Assembly::remove_unused_contigs(){
  cout << Date() << ": removing unused contigs and renumbering" << endl;
   vec<int> tigs_used( efastas.size(), 0);
  for ( size_t si = 0; si < scaffolds.size(); si++ ){
    const superb & s = scaffolds[si];
    ForceAssertGt( s.Ntigs(), 0 );
    for ( int tpos = 0; tpos < s.Ntigs(); tpos++ ){
      size_t tid = s.Tig(tpos);
      ForceAssertLt( tid, efastas.size() );
      if ( efastas[tid].Length1() != s.Len(tpos) ){
	PRINT5( si, tpos, tid, s.Len(tpos), efastas[tid].Length1() );
	ForceAssertEq( efastas[tid].Length1(), s.Len(tpos) );
      }
      tigs_used[tid]++;
    }
  }

  size_t unusedCt = 0;
  vec<int> offsets( efastas.size(), 0 );
  int offset = 0;
  for ( size_t tid = 0; tid < efastas.size(); tid++ ){
    if ( ! tigs_used[tid] ){
      offsets[tid] = -1;
      offset++;
      unusedCt++;
    }
    else{ offsets[tid] = offset; }
  }
  cout << Date( ) << ": found " << unusedCt << " unused contigs, removing" << endl;
  ForceAssertEq( offsets.size(), efastas.size() );
  for ( int tid = 0; tid < offsets.isize(); tid++ ){
    if ( offsets[tid] > 0 ){
      efastas[ tid - offsets[tid] ] = efastas[tid];
      tigMap[tid - offsets[tid] ] = tigMap[tid];
    }      
  }
  efastas.resize( efastas.size() - offset );
  tigMap.resize( tigMap.size() - offset );
  ForceAssertEq( efastas.size(), tigMap.size() );
  
  cout << Date() << ": updating scaffolds tig ids" << endl;
  for ( size_t si = 0; si < scaffolds.size(); si++ ){
    for ( int tpos = 0; tpos < scaffolds[si].Ntigs(); tpos++ ){
      int tid = scaffolds[si].Tig(tpos);
      if ( offsets[tid] >= 0 ){
	int newtid = tid - offsets[tid];
	ForceAssertGe( newtid, 0 );
	scaffolds[si].SetTig( tpos, newtid );
      }
      else{
	scaffolds[si].RemoveTigByPos( tpos );
	tpos--;
      }
    }
  }
  cout << Date() << ": done with removing unused contigs" << endl;
}


void Assembly::reorder(){
  // Sorting scaffolds according to size and renumbering contigs according
  //   to sequential appearance in scaffolds
  cout << Date() << " sorting scaffolds" << endl;
  SortSync( scaffolds, scaffMap, scompare );
  renumber();
}



// renumber all the contigs sequentially according to the scaffold
void Assembly::renumber(){
  cout << Date() << " renumbering contigs for ordered scaffolds" << endl;
  vec<String> otigMap;  
  vec<superb> oscaffolds = scaffolds;
  vec<efasta> oefastas;
  int cTid = -1;
  for ( size_t is = 0; is < scaffolds.size(); is++ ){
    for ( int tpos = 0; tpos < scaffolds[is].Ntigs(); tpos++ ){
      cTid++;
      int oTid = scaffolds[is].Tig(tpos);
      oscaffolds.at(is).SetTig( tpos, cTid );
      oefastas.push_back( efastas.at(oTid) );
      otigMap.push_back( tigMap.at(oTid) );
    }
  }
  efastas.resize(0);
  tigMap.resize(0);
  size_t newScaffoldsTotLen = 0, newScaffoldsRedLen = 0;
  for ( size_t is = 0; is < scaffolds.size(); is++ ){
    newScaffoldsTotLen += scaffolds[is].FullLength();
    newScaffoldsRedLen += scaffolds[is].ReducedLength();
  }
  scaffolds = oscaffolds; 
  oscaffolds.resize(0);
  efastas = oefastas; 
  oefastas.resize(0);
  tigMap = otigMap;
}

void Assembly::dedup() {
  // Remove duplicate scaffolds.
  cout << Date() << " removing duplicate scaffolds: " << endl;
  cout << "initial number of scaffolds = " << scaffolds.size() << endl;
  int removed_count = 0;
  int removed_size = 0;
  vec<longlong> s2minRedLen( scaffolds.size(), 0 ), s2maxRedLen( scaffolds.size(), 0 );
  vec<longlong> s2minFullLen( scaffolds.size(), 0 ), s2maxFullLen( scaffolds.size(), 0 );
  for ( int si = 0; si < scaffolds.isize(); si++ ){
    //PRINT4( si, scaffolds[si].Ntigs(), scaffolds[si].FullLength(), scaffolds[si].ReducedLength() );
    for ( int i = 0; i < scaffolds[si].Ntigs(); i++ ){
      int ti = scaffolds[si].Tig(i);
      efasta &tigi = efastas[ ti ];
      s2minRedLen[si] += tigi.MinLength();
      s2maxRedLen[si] += tigi.MaxLength();
      s2minFullLen[si] += tigi.MinLength();
      s2maxFullLen[si] += tigi.MaxLength();
      if ( i < scaffolds[si].Ntigs() -1 ){
	s2minFullLen[si] += scaffolds[si].Gap(i) - 3 * scaffolds[si].Dev(i); 
	s2maxFullLen[si] += scaffolds[si].Gap(i) + 3 * scaffolds[si].Dev(i); 
      }
    }
  }

  for ( int si = 0; si < scaffolds.isize(); si++ ) {
    int Ni = scaffolds[si].Ntigs();
    for ( int sj = si + 1; sj < scaffolds.isize(); sj++ ) {
      int Nj = scaffolds[sj].Ntigs();
      if ( Nj != Ni ) 
	continue;
      if ( s2maxRedLen[si] < s2minRedLen[sj] )
	continue;
      if ( s2maxRedLen[sj] < s2minRedLen[si] )
	continue;
      if ( s2maxFullLen[si] < s2minFullLen[sj] )
	continue;
      if ( s2maxFullLen[sj] < s2minFullLen[si] )
	continue;
      
      // check if reverse complement
      cout << " passed initial tests\n";
      Bool isRc = True; //initial condition
      for ( int i = 0; i < Ni; i++ ){
	
	int ti = scaffolds[si].Tig(i);
	int tj = scaffolds[sj].Tig( Nj -1 -i );
	
	efasta &tigi = efastas[ ti ];
	efasta &tigj = efastas[ tj ];
	
	vec<basevector> v_ibases, v_jbases;
	tigi.ExpandTo(v_ibases);
	tigj.ExpandTo(v_jbases);

	Bool foundEqual = False;
	for ( size_t vi = 0; vi < v_ibases.size() && ! foundEqual; vi++ ){
	  for ( size_t vj = 0; vj < v_jbases.size(); vj++ ){
	    if ( v_jbases[vj].size() != v_ibases[vi].size() )
	      continue;
	    basevector jbases = v_jbases[vj];
	    jbases.ReverseComplement();
	    if ( jbases == v_ibases[vi] ){
	      foundEqual = True;
	      break;
	    }
	  }
	}
	if ( ! foundEqual ){
	  isRc = False;
	  break;
	}
      }
      // check if exact duplicate
      Bool isExactDup = True; //initial condition
      if ( ! isRc ){
	for ( int i = 0; i < Ni; i++ ){
	  
	  int ti = scaffolds[si].Tig(i);
	  int tj = scaffolds[sj].Tig(i);
	  
	  efasta &tigi = efastas[ ti ];
	  efasta &tigj = efastas[ tj ];
	  
	  vec<basevector> v_ibases, v_jbases;
	  tigi.ExpandTo(v_ibases);
	  tigj.ExpandTo(v_jbases);
	  
	  Bool foundEqual = False;
	  for ( size_t vi = 0; vi < v_ibases.size() && ! foundEqual; vi++ ){
	    for ( size_t vj = 0; vj < v_jbases.size(); vj++ ){
	      if ( v_jbases[vj].size() != v_ibases[vi].size() )
		continue;
	      if ( v_jbases[vj] == v_ibases[vi] ){
		foundEqual = True;
		break;
	      }
	    }
	  }
	  if ( ! foundEqual ){
	    isExactDup = False;
	    break;
	  }
	}
      }else{
	isExactDup = False;
      }
      if ( isRc || isExactDup ){
	String type = isRc ? "reverse complement" : "exact duplicate";
	
	cout << "scaffold " << sj << " is " << type << " of " << si
	     << " (length " << scaffolds[si].FullLength() << ")" << endl;
	
	scaffolds.erase( scaffolds.begin() + sj );
	scaffMap.erase( scaffMap.begin() + sj );
	s2minRedLen.erase( s2minRedLen.begin() + sj );
	s2maxRedLen.erase( s2maxRedLen.begin() + sj );
	s2minFullLen.erase( s2minFullLen.begin() + sj );
	s2maxFullLen.erase( s2maxFullLen.begin() + sj );
	sj--;
	removed_count++;
	removed_size += scaffolds[sj].ReducedLength();
      }
    }
  }
  cout << "removed " << removed_count << " duplicate scaffolds"
       << " (" << removed_size << " bases)" << endl;
  cout << "final number of scaffolds = " << scaffolds.size() << endl;
  remove_unused_contigs();
}


void Assembly::dedup2() {
  // Remove scaffolds that are possibly duplicate
  // The criteria are:
  // 1. two scaffolds have same number of contigs (n_contig)
  // 2. n_contig >= 2
  // 3. Each contig and gap much match 
  //    - Gaps are matched whan [ gap_size +/- 3 * std ] overlap
  //    - Contigs are matched when 
  //      - perfect efasta match for contig size < 50,000
  //      - 1/100 mismatch kmer rate for contig size >= 50,000 (!!!!! this is only meant to be temporary fix to the problem)
  // 4. Both rc and fw duplicates are checked

  int VERBOSITY = 0;
  const int EfastaMatchSize = 50 * 1000; 
  const int MaxDev  = 3;
  const double MaxMismatchRate = 0.01;
  cout << Date() << ": " << "Remove possible duplicate scaffolds" << endl;
  cout << Date() << ": " << "initial number of scaffolds = " << scaffolds.size() << endl;
  
  int removed_count = 0;
  int removed_size = 0;
  vec <Bool> to_remove ( scaffolds.size(), False );
  for ( int si = 0; si < scaffolds.isize(); si++ ) {
    if ( to_remove[si] ) continue;
    int Ni = scaffolds[si].Ntigs();
    //if ( Ni < 2 ) continue;
    for ( int sj = si + 1; sj < scaffolds.isize(); sj++ ) {
      if ( to_remove[sj] ) continue;
      int Nj = scaffolds[sj].Ntigs();
      if ( Nj != Ni ) continue;
      //VERBOSITY = ( si == 1 && sj == 6 ? 1: 0 );
      // check scaffolds duplicate in two passes: fw in pass=0, rc in pass=1
      for( int pass = 0; pass < 2; pass++ ) {
	if ( VERBOSITY >= 1 )
	cout << Date() << ": " << "pass= " << pass << endl;
	if ( VERBOSITY >= 1 )
	  cout << Date() << ": " << "check gap " << endl ;
	// compare the gap of two scaffolds
	bool gap_match = true;
	for ( int igap = 0; igap < Ni-1; ++igap ) {
	  int gap1 = scaffolds[si].Gap(igap);
	  int dev1 = scaffolds[si].Dev(igap);
	  int gap2 = scaffolds[sj].Gap( pass == 0 ? igap : Ni -2 - igap);
	  int dev2 = scaffolds[sj].Dev( pass == 0 ? igap : Ni -2 - igap);
	  if ( IntervalOverlap( gap1 - MaxDev * dev1, gap1 + MaxDev * dev1 + 1,
		gap2 - MaxDev * dev2, gap2 + MaxDev * dev2 + 1) == 0 ) {
	    gap_match = false;
	    break;
	  }
	}
	if ( !gap_match) continue;
	if ( VERBOSITY >= 1 )
	  cout << Date() << ": " << "check contigs " << endl ;

	// check if the contigs matchs
	bool tig_match = true;
	if ( VERBOSITY >= 1 )
	  cout << "ncontigs= " << Ni << endl;
	for ( int itig = 0; itig < Ni; ++itig) {
	  if ( VERBOSITY >= 1 )
	    cout << " check contig " << itig << ": "
	    << " vs " << scaffolds[sj].Tig( pass == 0 ? itig :  Nj -1 -itig ) << endl;
	  efasta &tigi = efastas[ scaffolds[si].Tig(itig) ];
	  efasta &tigj = efastas[ scaffolds[sj].Tig( pass == 0 ? itig :  Nj -1 -itig ) ] ;

	  if ( VERBOSITY >= 1 )
	    cout << "tigi.size()= " << tigi.size() << endl;
	  if ( VERBOSITY >= 1 )
	    cout << "tigj.size()= " << tigj.size() << endl;
	  
	  // check if contig sizes match
	  int tig_size = ( tigi.size() + tigj.size()  ) /2;
	  int MaxMismatch =  tig_size * MaxMismatchRate ;
	  if ( abs( (int)tigi.size() - (int)tigj.size() ) > MaxMismatch ) { tig_match = false; break; }

	  // require perfect efasta match if contig size less than EfastaMatchSize
	  if ( tig_size < EfastaMatchSize ) {
	    if ( VERBOSITY >= 1 )
	      cout << " check efasta " << endl;
	    vec<basevector> v_ibases, v_jbases;
	    tigi.ExpandTo(v_ibases);
	    tigj.ExpandTo(v_jbases);
	    if ( pass == 1 ) 
	      for ( size_t k = 0; k < v_jbases.size(); ++k ) { v_jbases[k].ReverseComplement(); }
	    bool foundEqual = False;
	    for ( size_t vi = 0; vi < v_ibases.size() && ! foundEqual; vi++ ){
	      for ( size_t vj = 0; vj < v_jbases.size(); vj++ ){
		if ( v_jbases[vj].size() != v_ibases[vi].size() )
		  continue;
		basevector jbases = v_jbases[vj];
		if ( jbases == v_ibases[vi] ){
		  foundEqual = True;
		  break;
		}
	      }
	    }
	    if ( ! foundEqual ) {
	      tig_match = false;
	      break;
	    }
	  } 
	  // larger contig size. do kmer matching
	  else {
	    if ( VERBOSITY >= 1 )
	      cout << " check kmers " << endl;
	    basevector base1, base2;
	    tigi.FlattenTo( base1 );
	    tigj.FlattenTo( base2 );
	    if ( pass == 1 ) base2.ReverseComplement();
	    const int K = 24;
	    ForceAssertGt( base1.isize( ), K );
	    vec< basevector > kmers1( base1.isize( ) - K + 1);
            #pragma omp parallel for
	    for ( int jz = 0; jz <= base1.isize( ) - K; jz += 1000 ) 
	      for ( int j = jz; j <= Min( base1.isize( ) - K, jz + 1000 ); j++ ) 
		kmers1[j].SetToSubOf( base1, j, K ); 
	      ParallelUniqueSort(kmers1);    
	    ForceAssertGt( base2.isize( ), K );
	    vec< basevector > kmers2( base2.isize( ) - K + 1);
            #pragma omp parallel for
	    for ( int jz = 0; jz <= base2.isize( ) - K; jz += 1000 ) 
	      for ( int j = jz; j <= Min( base2.isize( ) - K, jz + 1000 ); j++ ) 
	        kmers2[j].SetToSubOf( base2, j, K ); 
	    ParallelUniqueSort(kmers2);    

	    // compare how many kmers are identical for the two sorted list
	    int nkmer1 = kmers1.size(), nkmer2 = kmers2.size();
	    int nkmer = (nkmer1 + nkmer2)/2;
	    int count = 0;
	    for( size_t i = 0, j = 0; i < kmers1.size() && j < kmers2.size(); ) {
	      if ( kmers1[i] > kmers2[j] ) j++;
	      else if ( kmers1[i] < kmers2[j] )  i++;
	      else count++, i++, j++;
	    }
	    if ( VERBOSITY >= 1 ) {
	      cout << "nkmer= " << nkmer << endl;
	      cout << "duplicate= " << count << endl;
	    }
	    if ( abs(nkmer - count) > int( nkmer * MaxMismatchRate) ) {
	      tig_match = false;
	      break;
	    }
	  }
	} // for itig
	if ( !tig_match ) continue;
	// ---------------------------------------------------------------------------
	// Now we concluded that the two scaffolds are duplicate. Remove the later one
	// --------------------------------------------------------------------------
	{
	  String type = pass == 0 ? "fw duplicate" : "rc duplicate";
	  cout << "scaffold " << sj << " is " << type << " of " << si
	    << " (length " << scaffolds[si].FullLength() << ")" << endl;
	  to_remove[sj] = True;
	  removed_count++;
	  removed_size += scaffolds[sj].ReducedLength();
	  break; // do not go second pass 
	}
      } // end pass 2
    } // end for sj
  } // end for si
  
  // now remove the duplicate scaffolds 
  {
    EraseIf( scaffolds, to_remove );
    EraseIf( scaffMap, to_remove );
  }

  cout << "removed " << removed_count << " duplicate scaffolds"
       << " (" << removed_size << " bases)" << endl;
  cout << "final number of scaffolds = " << scaffolds.size() << endl;
  remove_unused_contigs();
}


void Assembly::Write( const String head_out ) const {

  // writing output
  cout << Date() << ": writing output files" << endl;
  WriteSuperbs( head_out + ".superb", scaffolds );
  WriteSummary( head_out + ".summary", scaffolds );


  
  Ofstream( efout, head_out + ".contigs.efasta" );
  for ( size_t id = 0; id < efastas.size(); id++ )
    efastas[id].Print(efout, "contig_" + ToString(id) );
}

void Assembly::WriteExtra( const String head_out ) const{

  vec<fastavector> fastas(efastas.size());
  for ( size_t id = 0; id < efastas.size(); id++ )
    efastas[id].FlattenTo( fastas[id] );
  Ofstream( fout, head_out + ".contigs.fasta" );
  for ( size_t id = 0; id < fastas.size(); id++ )
    fastas[id].Print(fout, "contig_" + ToString(id) );

  {
    vecfastavector vec_tigs;
    for ( size_t i = 0; i < fastas.size( ); i++ )
      vec_tigs.push_back_reserve( fastas[i] );
    vec_tigs.WriteAll( head_out + ".contigs.vecfasta" ); 
  }

  vecbasevector bases( efastas.size() );
  for ( size_t id = 0; id < efastas.size(); id++ )
    efastas[id].FlattenTo( bases[id] );
  
  bases.WriteAll( head_out + ".contigs.fastb" );
  
  Ofstream( cmout, head_out + ".contigs.mapping" );
  for ( size_t id = 0; id < tigMap.size(); id++ )
    cmout << ToString(id) + " from " + ToString( tigMap[id] ) << "\n";

  Ofstream( smout, head_out + ".superb.mapping" );
  for ( size_t is = 0; is < scaffMap.size(); is++ )
    smout << ToString(is) + " from " + ToString( scaffMap[is] ) << "\n";

  WriteScaffoldedEFasta( head_out + ".assembly.efasta", efastas, scaffolds );
  WriteScaffoldedFasta( head_out + ".assembly.fasta", fastas, scaffolds );

}
