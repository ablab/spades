///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2012) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

const char *DOC =
  "Prepares assembly for submission to NCBI (or similar). It removes small "
  "contigs and scaffolds, removes duplicate scaffolds, sets a minimum gap size, "
  "generates AGP and tbl files, and optionally removes identified contaminants.";


#include "CoreTools.h"
#include "MainTools.h"
#include "TokenizeString.h"
#include "agp/AgpFile.h"
#include "efasta/EfastaTools.h"
#include "VecUtilities.h"
#include "paths/AssemblyCleanupTools.h"

/**
 * SubmissionPrep
 *
 * Prepare scaffolds for submission (order, clean, remove small contigs, etc. ).
 *
 */

size_t scaffoldsTotLen( const vec<superb>& scaffolds ){
  size_t len = 0;
  for ( size_t is = 0; is < scaffolds.size(); is++ )
    len += scaffolds[is].FullLength();
  return len;
}

size_t scaffoldsRedLen( const vec<superb>& scaffolds ){
  size_t len = 0;
  for ( size_t is = 0; is < scaffolds.size(); is++ )
    len += scaffolds[is].ReducedLength();
  return len;
}


String zero_padded_int(const int value, const int length) {
  ostringstream out;
  out.fill('0');
  out.width(length);
  out << value;
  return out.str();
}

String contig_id_str(const int contig_id)
{
  return "contig" + zero_padded_int(contig_id+1, 6);
}

// restrict efasta given start position and length measured
// using the first option in the curly brackets.
void RestrictEfasta( efasta& source, size_t start, size_t len ){
  int nstart = source.Index1(start);
  int nend = source.Index1(start+len);

  ForceAssert( nend >= nstart );
  ForceAssert( nend <= source.isize() );
  
  source.erase( nend, source.size() - nend );
  source.erase( 0, nstart );
  return;
}

  
void set_min_gaps(vec<superb>& scaffolds, const int min_gap)
{
  cout << Date() << " resetting gaps < " << min_gap << endl;
  for ( size_t is = 0; is < scaffolds.size(); is++ )
    for ( int tpos = 0; tpos < scaffolds[is].Ntigs() -1; tpos++ )
      if ( scaffolds[is].Gap(tpos) < min_gap ) {
	scaffolds[is].SetGap( tpos, min_gap );
	scaffolds[is].SetDev( tpos, 0 );
      }
}


int load_contamination_list( String in_contam_file, vec<vec<triple<int,int,int> > >& v_contam, const vec<efasta>& efastas ){
  if ( ! IsRegularFile( in_contam_file ) )
      FatalErr("contamination file " + in_contam_file + " not found");
  ifstream in( in_contam_file.c_str() );
  String line;
  int contamTotLen = 0;
  // reading data in the format:
  cout << Date() << " reading in contamination file" << flush;
  while ( getline(in,line) ){
    vec<String> tokens;
    Tokenize( line, tokens );
    // Ignore anything beyond 4th column
    ForceAssert( tokens.size() >= 4);
    int tid    = tokens[0].Int();
    int tlen   = tokens[1].Int();
    int cbegin = tokens[2].Int();
    int cend   = tokens[3].Int();
    ForceAssertGe( tid, 0 );
    ForceAssertLt( tid, efastas.isize() );
    ForceAssertEq( tlen, efastas[tid].Length1() );
    ForceAssertGe( cbegin, 0 );
    ForceAssertLe( cend, efastas[tid].Length1() );
    contamTotLen += cend - cbegin;
    v_contam.at(tid).push_back( triple<int,int,int>(tlen,cbegin, cend) );
  }
  in.close();
  cout << " (length = " << contamTotLen << ")" << endl;
  return contamTotLen;
}

int remove_contamination_list( const vec< vec<triple<int,int,int> > >& v_contam,
			       vec<efasta>& efastas, vec<superb>& scaffolds,
			       vec<String>& tigMap, const String save_contam_file ){


  vec<superb> new_tscaffolds( efastas.size() ); 
  vec<Bool> modif_contigs( efastas.size(), False );
  vec< pair<int,int> > beg_gaps( efastas.size() );
  vec< pair<int,int> > end_gaps( efastas.size() );
  vec<fastavector> contamination;
  vec<String> contamination_ids;

  int contamCheckLen = 0;
  for ( size_t tid = 0; tid < v_contam.size(); tid++ ){  
    if ( v_contam.at(tid).size() == 0 ) 
      continue;
    else
      modif_contigs[tid] = True;
    //cout << "\n  -------- \n\n"; PRINT(tid);

    basevector tbases;
    efastas[tid].FlattenTo( tbases );

    int specContamLen = 0, resContamLen = 0;
    vec<Bool> treg( efastas[tid].Length1() +2, True );
    for ( size_t i = 0; i < v_contam[tid].size(); i++ ){
      int clen = v_contam[tid][i].first;
      int cbeg = v_contam[tid][i].second;
      int cend = v_contam[tid][i].third;
      ForceAssertEq( (int)efastas[tid].Length1(), clen );
      specContamLen +=  cend - cbeg;
      for ( int j = cbeg; j < cend; j++ )
	treg[j+1] = False;

      basevector contam_bv(tbases, cbeg, specContamLen);
      fastavector contam_fv(contam_bv);
      contamination.push_back(contam_fv);
      contamination_ids.push_back("contam_contig" + zero_padded_int(tid, 6) + "_" + ToString(cbeg) + "_" + ToString(cend));
    }
    for ( size_t i = 1; i < treg.size() -1; i++ )
      if ( ! treg[i] ) {
	resContamLen++;
	contamCheckLen++;
      }
    ForceAssertEq( specContamLen, resContamLen );
    
    treg.front() = treg.back() = False; 
    
    vec<int> begs, ends;
    
    for ( size_t p = 1; p < treg.size(); p++ )
      if ( ! treg[p -1]  && treg[p] )
	begs.push_back( p - 1 );
      else if ( treg[p -1] && ! treg[p] )
	ends.push_back( p - 1 );
    
    ForceAssert( begs.size() == ends.size() );
    if ( begs.size() > 0 ){
      if ( begs.front() > 0 ){
	beg_gaps[tid].first  = begs[0];
	beg_gaps[tid].second = 1;
      }
      if ( ends.back() != efastas[tid].Length1() ){
	end_gaps[tid].first  = efastas[tid].Length1() - ends.back();
	end_gaps[tid].second = 1;
      }
    }else{
      // entire contig removed      
      int b = efastas[tid].Length1() / 2;
      int e = b + efastas[tid].Length1()%2;
      beg_gaps[tid] = pair<int,int>( b, 1 );
      end_gaps[tid] = pair<int,int>( e, 1 );
    }
    
    //cout << "updating sequences" << endl;
    
    efasta tefasta    = efastas[tid];
    new_tscaffolds[tid].SetNtigs( begs.size() );
    if ( begs.size() == 0 ){
      efastas[tid].resize(0);
      new_tscaffolds[tid].SetNtigs( 1 );
      new_tscaffolds[tid].SetTig(0, tid);
      new_tscaffolds[tid].SetLen(0, 0);
    }else {
      for ( size_t i = 0; i < begs.size(); i++ ){
	basevector lbases( tbases, begs[i], ends[i] - begs[i] );
	efasta lefasta = tefasta;
	//cout << "restricting efasta" << endl;
	RestrictEfasta( lefasta, begs[i], ends[i] - begs[i] );
	basevector tmpbases;
	lefasta.FlattenTo(tmpbases);
	ForceAssert( lbases == tmpbases );
	new_tscaffolds[tid].SetLen(i, lbases.size());
	//cout << "updating efasta and base vecs" << endl;
	if ( i == 0 ){
	  new_tscaffolds[tid].SetTig(i, tid);
	  efastas[tid] = lefasta;
	}else{
	  new_tscaffolds[tid].SetTig(i, efastas.size() );
	  efastas.push_back( lefasta );
	  tigMap.push_back( ToString(tid) );
	}
	if ( i != begs.size() -1 ){
	  new_tscaffolds[tid].SetGap(i, begs[i+1] - ends[i]);
	  new_tscaffolds[tid].SetDev(i, 1);
	}
      }
    }
  }
  
  // update scaffolds
  size_t endGapSum = 0;
  for ( size_t si = 0; si < scaffolds.size(); si++ ){
    superb & s = scaffolds[si];
    for ( int tpos = 0; tpos < s.Ntigs(); tpos++ ){
      size_t tid = s.Tig(tpos);
      if ( tid < modif_contigs.size() && modif_contigs[tid] ){
	//cout << "\n ------------\n";
	//new_tscaffolds[tid].PrintRange(cout, "insert:", 0, new_tscaffolds[tid].Ntigs() );
	size_t lenBefore = s.FullLength();
	size_t lenAfter = 0;
	if ( tpos == 0 ){
	  lenAfter += beg_gaps[tid].first;
	  if ( s.Ntigs() > 1 && new_tscaffolds[tid].FullLength() == 0
	       && new_tscaffolds[tid].Ntigs() == 0 ){
	    lenAfter += end_gaps[tid].first;
	    lenAfter += s.Gap(tpos);
	  }
	}
	if ( tpos == s.Ntigs() -1 ) {
	  lenAfter += end_gaps[tid].first;
	  if ( s.Ntigs() > 1 && new_tscaffolds[tid].FullLength() == 0
	       && new_tscaffolds[tid].Ntigs() == 0 ){
	    lenAfter += beg_gaps[tid].first;
	    lenAfter += s.Gap(tpos -1 );
	  }
	}
	endGapSum += lenAfter;
	//s.PrintRange( cout, "before:", tpos -1, tpos +1 );
	s.ReplaceTigBySuper( tpos, new_tscaffolds[tid], 
			     beg_gaps[tid].first, beg_gaps[tid].second,
			     end_gaps[tid].first, end_gaps[tid].second);
	//cout << "\nAfter" << endl;
	modif_contigs[tid] = False;
	//s.PrintRange( cout, "after:", tpos -1, tpos + new_tscaffolds[tid].Ntigs() );
	lenAfter += s.FullLength();
	ForceAssertEq( lenBefore, lenAfter );
	tpos += new_tscaffolds[tid].Ntigs() -1;
      }
    }
  }  
  if (save_contam_file != "") {
    cout << Date() << " saving contamination segments to " << save_contam_file << endl;
    Ofstream(contam_stream, save_contam_file);
    for (u_int i = 0; i < contamination.size(); ++i)
      contamination[i].Print(contam_stream, contamination_ids[i]);
    contam_stream.close();
  }
  return contamCheckLen;
}

void
translate_rings_file(const String input_filename, const String output_filename)
{
  Ifstream (in, input_filename);
  Ofstream (out, output_filename);
  String line;
  out << "#scaffold\tgap_size\tgap_dev\tn_links" << endl;
  while (getline(in,line)){
    vec<String> tokens;
    Tokenize( line, tokens );
    if (tokens.size() == 4 && tokens[0][0] >= '0' && tokens[0][0] <= '9') {
      int scaffold    = tokens[0].Int();
      int gap_size    = tokens[1].Int();
      int gap_dev     = tokens[2].Int();
      int n_links     = tokens[3].Int();

      out << "scaffold";
      out.width(5);
      out.fill('0');
      out << scaffold+1;
      out << "\t" << gap_size << "\t" << gap_dev << "\t" << n_links << endl;      
    }
  }
  in.close();
  out.close();
}


void
remove_contamination(String in_contam_file,
		     Assembly& assembly,
		     int min_contig_solo,
		     int min_contig_in,
		     String out_save_contam_file)
{
  vec<String>& scaffMap = assembly.scaffMap;
  vec<String>& tigMap = assembly.tigMap;
  vec<superb>& scaffolds = assembly.scaffolds;
  vec<efasta>& efastas = assembly.efastas;

  cout << Date() << " loading contaminant information" << endl;
  vec< vec<triple<int,int,int> > > v_contam( efastas.size() );
  int contamTotLen = 
    load_contamination_list( in_contam_file, v_contam, efastas );
  int contamCheckLen = 
    remove_contamination_list( v_contam, efastas, scaffolds, tigMap, out_save_contam_file);
  if ( contamTotLen != contamCheckLen ){
    cout << " There seems to be a problem with contamination specification, are there contaminant overlaps?" << endl;
    ForceAssertEq( contamTotLen, contamCheckLen );
  }
  assembly.remove_small_contigs(min_contig_solo, min_contig_in);
  assembly.remove_small_scaffolds(min_contig_solo);
  assembly.remove_unused_contigs();
  assembly.check_integrity();
}


int main( int argc, char *argv[] )
{
  RunTime( );
  
  BeginCommandArguments;
  CommandDoc(DOC);
  CommandArgument_String_Doc( HEAD_IN,
    "Input assembly files <HEAD_IN>.*");
  CommandArgument_String_Doc( HEAD_OUT,
    "Output assembly files <HEAD_OUT>.*" );
  CommandArgument_Int_OrDefault_Doc( MIN_CONTIG_SIZE_SOLO, 1000, 
    "Remove solo contigs that are smaller than MIN_CONTIG_SIZE_SOLO");
  CommandArgument_Int_OrDefault_Doc( MIN_CONTIG_SIZE_IN, 200, 
    "Remove contigs within scaffolds that are smaller than MIN_CONTIG_SIZE_IN");
  CommandArgument_Bool_OrDefault_Doc( REORDER, False,
    "Reorder scaffolds, largest to smallest.");
  CommandArgument_Bool_OrDefault_Doc( DEDUP, True,
    "Remove duplicate scaffolds.");
  CommandArgument_Int_OrDefault_Doc( MIN_GAP, 100, 
    "Gaps smaller than this will be increased to this value (per NCBI)");
  CommandArgument_String_OrDefault_Doc( CONTAM, "", 
    "Full path name of contamination file with: contig_id contig_lenght begin end");
  CommandArgument_Bool_OrDefault_Doc( REMOVE_CONTAM_LAST, True,
    "Remove contamination information at the end of the process. This option assumes that contamination coordinates corresopond to the modified assembly created by the first part of the code.");
  CommandArgument_Bool_OrDefault_Doc( SAVE_CONTAM, True,
    "Save fasta of removed contamination segments.");
  CommandArgument_UnsignedInt_OrDefault_Doc( FASTA_BATCH_SIZE, 10000,
    "Output contig fastas won't contain more than this many elements.");
  CommandArgument_Bool_OrDefault_Doc( AMBIGUOUS_BASE_CODES, False,
    "Allow IUPAC ambiguous base codes in contigs.  If false, SNP variations will be annotated in tbl file.");
  CommandArgument_String_OrDefault_Doc( ISOLATE_CONTIGS, "",
    "Detach specified contigs from their scaffolds, but don't remove them.");
  CommandArgument_String_OrDefault_Doc( EXTRA_CONTIGS, "",
    "Fasta file with extra contigs to be merged into the assembly as singleton contigs.");
  EndCommandArguments;

  // File names.
  String in_efasta_file = HEAD_IN + ".contigs.efasta";
  String in_superb_file = HEAD_IN + ".superb";
  String in_contam_file = CONTAM;

  String out_efasta_file        = HEAD_OUT + ".contigs.efasta";
  String out_fasta_file         = HEAD_OUT + ".contigs.fasta";
  String out_tbl_file           = HEAD_OUT + ".contigs.tbl";
  String out_superb_file        = HEAD_OUT + ".superb";
  String out_fastb_file         = HEAD_OUT + ".fastb";
  String out_contig_map_file    = HEAD_OUT + ".contigs.mapping";
  String out_scaffold_map_file  = HEAD_OUT + ".superb.mapping";
  String out_save_contam_file   = SAVE_CONTAM ? HEAD_OUT + ".contam.fasta" : "";

  vec<int> isolate_contigs;
  ParseIntSet( ISOLATE_CONTIGS, isolate_contigs, false );

  // Lodading scaffolds
  cout << Date( ) << " loading superb file" << endl;
  vec<superb> scaffolds;
  ReadSuperbs( in_superb_file, scaffolds );
  
  // reading contig information
  cout << Date( ) << " loading contigs efasta file" << endl;
  if ( ! IsRegularFile( in_efasta_file) )
    FatalErr("input file " + in_efasta_file + " not found");
  vec<efasta> efastas;
  LoadEfastaIntoStrings(in_efasta_file, efastas);

  if (EXTRA_CONTIGS != "") {
    cout << Date( ) << " loading extra contigs" << endl;
    vec<fastavector> extra_contigs;
    LoadFromFastaFile( EXTRA_CONTIGS, extra_contigs );
    
    // Add the extra contigs in as supers
    for (size_t i = 0; i < extra_contigs.size(); ++i) {
      size_t tig = efastas.size();
      efasta contig(extra_contigs[i]);
      superb super;
      super.PlaceFirstTig(tig, contig.Length1());
      scaffolds.push_back(super);
      efastas.push_back(contig);
    }
    cout << Date( ) << " " << extra_contigs.size() << " extra contigs appended as singleton scaffolds" << endl;
  }

  // make sure that there are no N's at the ends
  for ( size_t id = 0; id < efastas.size(); id++ ){
    basevector bases;
    efastas[id].FlattenTo( bases );
    String sbases = bases.ToString();
    ForceAssert( sbases.at(0) != 'N' );
    ForceAssert( sbases.at(0) != 'n' );

    ForceAssert( sbases.back() != 'N' );
    ForceAssert( sbases.back() != 'n' );
  }

  if (isolate_contigs.size() > 0) {
    int nscaffolds_orig_size = scaffolds.size();
    for (int i = 0; i < isolate_contigs.isize(); ++i) {
      int isolate = isolate_contigs[i];
      cout << Date( ) << " isolating contig " << isolate << " into it's own scaffold" << endl;
      for (int s = 0; s < nscaffolds_orig_size; ++s) {
	for (int t = 0; t < scaffolds[s].Ntigs(); ++t) {
	  int tig = scaffolds[s].Tig(t);
	  if (tig == isolate) {
	    int len = scaffolds[s].Len(t);
	    PRINT4(tig, s, t, len);
	    superb iso_super;
	    iso_super.PlaceFirstTig(tig, len);
	    scaffolds[s].RemoveTigByPos(t);
	    scaffolds.push_back(iso_super);
	    break;
	  }
	}
      }
    }
  }

  size_t origScaffoldsTotLen = scaffoldsTotLen( scaffolds );
  size_t origScaffoldsRedLen = scaffoldsRedLen( scaffolds );

  Assembly assembly( scaffolds, efastas );

  // check initial basic integrity
  assembly.check_integrity();

  vec<String> & scaffMap = assembly.scaffMap;
  vec<String> & tigMap = assembly.tigMap;

  // remove contamination
  vec<fastavector> contamination;
  if ( in_contam_file.nonempty() && !REMOVE_CONTAM_LAST) {
    remove_contamination(in_contam_file, assembly, MIN_CONTIG_SIZE_SOLO, MIN_CONTIG_SIZE_IN,
			 out_save_contam_file);
  }

  assembly.remove_small_contigs( MIN_CONTIG_SIZE_SOLO, MIN_CONTIG_SIZE_IN );
  assembly.remove_small_scaffolds( MIN_CONTIG_SIZE_SOLO );
  assembly.remove_unused_contigs();
  assembly.check_integrity();


  if (DEDUP) assembly.dedup();
  if (REORDER) assembly.reorder();

  // check integrity
  assembly.check_integrity();

  // remove contamination
  if ( in_contam_file.nonempty() && REMOVE_CONTAM_LAST) {
    remove_contamination(in_contam_file, assembly, MIN_CONTIG_SIZE_SOLO, MIN_CONTIG_SIZE_IN,
			 out_save_contam_file);
  }

  // reset small or negative gaps
  set_min_gaps(assembly.scaffolds, MIN_GAP);
  

  // writing output
  cout << Date() << " writing superb" << endl;
  WriteSuperbs( out_superb_file, assembly.scaffolds );

  cout << Date() << " writing assembly efasta" << endl;
  Ofstream( efout, out_efasta_file );
  Ofstream( fout, out_fasta_file );
  Ofstream( tblout, out_tbl_file );

  cout << Date() << " writing fsa & tbl files" << endl;
  int batch = 0;
  String batch_suffix = (assembly.efastas.size()>FASTA_BATCH_SIZE) ? ("_" + zero_padded_int(batch, 2)) : "";
  Ofstream(fsa, HEAD_OUT + batch_suffix + ".fsa");
  Ofstream(tbl, HEAD_OUT + batch_suffix + ".tbl");

  vec<fastavector> fastas(assembly.efastas.size());
  int snp_count = 0,indel_count = 0, ambiguous_bases = 0;

  for ( size_t id = 0; id < assembly.efastas.size(); id++ ) {
    String name = contig_id_str(id);

    // output efasta record
    assembly.efastas[id].Print(efout, name);

    // chunk into files of no more than 10,000 contigs for NCBI
    if ((id+1) % FASTA_BATCH_SIZE == 0) {
      cout << "." << flush;
      fsa.close();
      ++batch;
      batch_suffix = "_" + zero_padded_int(batch, 2);
      OpenOfstream( fsa, HEAD_OUT + batch_suffix + ".fsa");
      tbl.close();
      OpenOfstream( tbl, HEAD_OUT + batch_suffix + ".tbl");
    }

    // code fragment from EfastaToFasta. --bruce
    // flatten to fsa (fasta)
    fastavector v;
    vec<Ambiguity> va;
    assembly.efastas[id].FlattenTo( v, va, AMBIGUOUS_BASE_CODES );
    fastas[id] = v;

    v.Print(fsa, name);
    v.Print(fout, name);

    if (va.size()) {
      ostringstream feature;

      feature << ">Feature " << name << "\n";

      // add annotation of alternate paths to tbl file
      for (size_t i = 0; i < va.size(); i++)
	feature << va[i].to_annotation();
      tbl << feature.str();
      tblout << feature.str();
    }

    int snps = 0, indels = 0;
    int amb =  assembly.efastas[id].AmbCount(snps, indels);
    ambiguous_bases += amb;
    snp_count += snps;
    indel_count += indels;
  }
  efout.close();
  tbl.close();
  tblout.close();
  fsa.close();
  fout.close();

  if (assembly.efastas.size() >= FASTA_BATCH_SIZE)
    cout << endl;

  // Generate agp file...code fragment from agp/SuperbsToAgp. --bruce
  cout << Date() << " writing agp file" << endl;
  Ofstream(agpfile, HEAD_OUT + ".agp");
  for (u_int ii=0; ii<assembly.scaffolds.size( ); ii++) {
    agp_chromosome agp( zero_padded_int(ii+1, 5 )) ;
    
    // Loop over the contigs.
    for (int jj=0; jj<assembly.scaffolds[ii].Ntigs( ); jj++) {

      // Add contig.
      String cg_name = contig_id_str(assembly.scaffolds[ii].Tig(jj));

      int cg_len = assembly.scaffolds[ii].Len( jj );
      ForceAssertEq((u_int)cg_len, fastas[assembly.scaffolds[ii].Tig( jj )].size());
      int cg_start = 0;
      int cg_stop = cg_len - 1;
      bool cg_RC = false;
      agp_contig new_contig( cg_name, cg_len, cg_start, cg_stop, cg_RC );
      new_contig.SetType( agp_contig::wgs_contig );
      agp.AddContig ( new_contig );

      // Add gap.
      if ( jj < assembly.scaffolds[ii].Ntigs( ) - 1 ) {
	String gap_type = "fragment";
	bool is_bridged = true;
	int gap_len = assembly.scaffolds[ii].Gap( jj );
	agp_gap new_gap( gap_type, gap_len, is_bridged );
	agp.AddGap( new_gap );
      }
    }
    
    // Print info for super.
    String base_name = "scaffold";
    agp.Print( agpfile, &base_name );
  }

  //vecbasevector obases( assembly.efastas.size() );
  //for ( size_t id = 0; id < assembly.efastas.size(); id++ )
  //  assembly.efastas[id].FlattenTo( obases[id] );
  
  //obases.WriteAll( out_fastb_file );


  cout << Date() << " writing gapped assembly efasta" << endl;
  WriteScaffoldedEFasta( HEAD_OUT + ".assembly.efasta", assembly.efastas, assembly.scaffolds, True);

  cout << Date() << " writing gapped assembly fasta" << endl;
  vec<Bool> rc(fastas.size(), False);
  WriteScaffoldedFasta( HEAD_OUT + ".assembly.fasta", fastas, assembly.scaffolds, rc, MIN_GAP, 'N', True);

  Ofstream( cmout, out_contig_map_file );
  for ( size_t id = 0; id < assembly.tigMap.size(); id++ )
    cmout << ToString(id) + " from " + ToString( assembly.tigMap[id] ) << "\n";

  Ofstream( smout, out_scaffold_map_file );
  for ( size_t is = 0; is < assembly.scaffMap.size(); is++ )
    smout << ToString(is) + " from " + ToString( assembly.scaffMap[is] ) << "\n";

  {
    size_t newScaffoldsTotLen = 0, newScaffoldsRedLen = 0, newContigTotLen = 0;
    vec<int> scaffolds_full(assembly.scaffolds.size());
    vec<int> scaffolds_reduced(assembly.scaffolds.size());
    for ( size_t is = 0; is < assembly.scaffolds.size(); is++ ){
      int tlen = assembly.scaffolds[is].FullLength();
      int rlen = assembly.scaffolds[is].ReducedLength();
      scaffolds_full[is] = tlen;
      newScaffoldsTotLen += tlen;
      scaffolds_reduced[is] = rlen;
      newScaffoldsRedLen += rlen;
    }
    vec<int> contig_sizes(fastas.size());
    for ( size_t ic = 0; ic < fastas.size(); ic++ ){
      int len = fastas[ic].size();
      newContigTotLen += len;
      contig_sizes[ic] = len;
    }

    String rings_file_in = HEAD_IN + ".rings";
    String rings_file_out = HEAD_OUT + ".rings";
    if (IsRegularFile(rings_file_in)) {
      cout << Date() << " Translating rings file" << endl;
      translate_rings_file(rings_file_in, rings_file_out);
    }


    cout << Date() << " Assembly statistics:" << endl;
    cout << "Number of scaffolds: " << ToStringAddCommas(assembly.scaffolds.size()) << endl;
    cout << "Total scaffold length including gaps: " << ToStringAddCommas(newScaffoldsTotLen)
	 << " (N50 = " << ToStringAddCommas(N50(scaffolds_full)) << ")" << endl;
    cout << "Total scaffold length excluding gaps: " << ToStringAddCommas(newScaffoldsRedLen)
	 << " (N50 = " << ToStringAddCommas(N50(scaffolds_reduced)) << ")" << endl;
    cout << "Number of contigs: " << ToStringAddCommas(fastas.size())
	 << " (N50 = " << ToStringAddCommas(N50(contig_sizes)) << ")" << endl;

    if (ambiguous_bases > 0)
      cout << "Ambiguous bases: " << ToStringAddCommas(ambiguous_bases)
	   << " (rate = 1/" << ToStringAddCommas(newScaffoldsRedLen / ambiguous_bases) << ")" << endl;
    if (snp_count > 0)
      cout << "SNP count: " << ToStringAddCommas(snp_count)
	   << " (rate = 1/" << ToStringAddCommas(newScaffoldsRedLen / snp_count) << ")" << endl;
    if (indel_count > 0)
      cout << "Indel count: " << ToStringAddCommas(indel_count)
	   << " (rate = 1/" << ToStringAddCommas(newScaffoldsRedLen / indel_count) << ")" << endl;
  }

  cout << Date() << " Done!" << endl;
}

