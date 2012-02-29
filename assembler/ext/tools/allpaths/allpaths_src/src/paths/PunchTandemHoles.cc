///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2011) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// PunchTandemHoles.  Find flagged, uncertain tandem repeats in an assembly.
// Remove them, preparatory to an attempt to reclose them with LongReadPostPatcher.
// This is to be followed by RecoverPunchTandem.

#include "Basevector.h"
#include "Fastavector.h"
#include "FastIfstream.h"
#include "MainTools.h"
#include "Superb.h"
#include "efasta/EfastaTools.h"
#include "paths/AssemblyCleanupTools.h"
#include "paths/AssemblyEdit.h"



class Edit
{
public:
  int i0;
  int i1;
  efasta efa;

  Edit(const int _i0, const int _i1, const efasta & _efa)
    : i0(_i0), i1(_i1), efa(_efa) 
  {
    ValidateEfastaRecord(_efa);
  }

  friend
  bool operator < (const Edit & a, const Edit & b) 
  { 
    if (a.i0 < b.i0) return true;
    if (a.i0 > b.i0) return false;
    if (a.i1 < b.i1) return true;
    if (a.i1 > b.i1) return false;
    return (a.efa < b.efa);
  }
};





void markups_read(const String & SCAFFOLDS_IN, 
                  vec<triple<int,int,int> > * markups_p,
                  const vec<efasta> & tigs_orig,
                  const int TIG,
                  const Bool VERBOSE)
{     
  if (IsRegularFile(SCAFFOLDS_IN + ".markup")) {
    cout << Date() << ": loading markups" << endl;
    fast_ifstream in(SCAFFOLDS_IN + ".markup");
    String line;
    while (1) {
      getline(in, line);
      if (in.fail()) break;
      int it; 
      int start;
      int stop;
      istringstream iline(line.c_str());
      iline >> it >> start >> stop;
      if (TIG < 0 || it == TIG) {
        if (VERBOSE) {
          cout << "markup: ";
          String s = tigs_orig[it].substr(start, stop - start);
          PRINT4(it, start, stop, s);
        }
        markups_p->push(it, start, stop);
      }
    }
  }
  cout << Date() << ": found " << markups_p->size() << " markups" << endl;
}






void edits_compile(const efasta & tig_orig, 
                   const vec< triple<int,int,int> > & markups,
                   const size_t it,
                   vec<Edit> * edits_p)
{
           
  // Make a list of the ambiguities in the efasta contig.
           
  for (int i = 0; i < tig_orig.isize(); i++) {    
    if (tig_orig[i] == '{') {    
      const int i0 = i;
      for (i++; i < tig_orig.isize(); i++) if (tig_orig[i] == ',') break;    
      for (i++; i < tig_orig.isize(); i++) if (tig_orig[i] == '}') break;
      const int i1 = i + 1;
      edits_p->push(i0, i1, tig_orig.substr(i0, i1 - i0)); 
    }   
  }
           
  // Add in markups.
           
  for (size_t im = 0; im < markups.size(); im++) {    
    if (markups[im].first == int(it)) {
      const int i0 = markups[im].second;
      const int i1 = markups[im].third;
      edits_p->push(i0, i1, tig_orig.substr(i0, i1 - i0));
    }    
  }
           
           
  // Sort
  Sort(*edits_p);
}



// Handle true tandems (not markups).

bool handle_true_tandem(const efasta & tig_orig,
                        Edit * edit_p,
                        const int MIN_PERIOD,
                        const bool USE_TANDEM)
{         
  if (edit_p->efa.Contains("{")) {
    vec<basevector> bvv;
    edit_p->efa.ExpandTo(bvv);
    sort(bvv.begin(), bvv.end(), LtBySize<basevector>());
    if (bvv[0].size() != 0) return false;
    String repeat;
    String bv1 = bvv[1].ToString();
    for (int div = bvv[1].isize(); div >= 1; div--) {
      if (bv1.isize() % div != 0) continue;
      int n = bv1.isize() / div;
      repeat = bv1.substr(0, n);
      String whole;
      for (int j = 0; j < div; j++)
        whole += repeat;
      if (whole == bv1) break;
    }
    bool tandem = true;
    for (int j = 2; j < bvv.isize(); j++) {
      String whole;
      for (int l = 0; l < bvv[j].isize() / repeat.isize(); l++)
        whole += repeat;
      if (bvv[j].ToString() != whole) tandem = false;
    }
    if (!tandem || repeat.isize() < MIN_PERIOD) return false;
    if (!USE_TANDEM) return false;
    while (tig_orig.Contains(repeat, edit_p->i0 - repeat.isize())) {
      edit_p->i0 -= repeat.isize();
      edit_p->efa = repeat + edit_p->efa;
    }
    while (tig_orig.Contains(repeat, edit_p->i1)) {
      edit_p->i1 += repeat.isize();
      edit_p->efa = edit_p->efa + repeat;
    }
  }
  return true;
}
         










int main(int argc, char *argv[])
{
  RunTime();

  BeginCommandArguments;
  CommandArgument_String(SCAFFOLDS_IN);
  CommandArgument_String_OrDefault(SCAFFOLDS_OUT, SCAFFOLDS_IN + ".tpunch");
  CommandArgument_Int_OrDefault(MIN_PERIOD, 5);
  CommandArgument_Int_OrDefault(FLANK, 100);
  CommandArgument_Bool_OrDefault(WRITE, True);
  CommandArgument_Bool_OrDefault(VERBOSE, False);
  CommandArgument_Int_OrDefault_Doc(TIG, -1, 
                                    "if specified, use only this contig");
  CommandArgument_Bool_OrDefault(USE_MARKUPS, True);
  CommandArgument_Bool_OrDefault(USE_TANDEM, True);
  EndCommandArguments;

  // Load assembly.

  cout << Date() << ": loading data" << endl;
  String supers_file = SCAFFOLDS_IN + ".superb";
  String efasta_file = SCAFFOLDS_IN + ".contigs.efasta";
  vec<superb> scaffolds;
  ReadSuperbs(supers_file, scaffolds);
  vec<efasta> tigs_orig;
  LoadEfastaIntoStrings(efasta_file, tigs_orig);

  const size_t nt_orig = tigs_orig.size();

     
  // validate contigs

  for (size_t it = 0; it < nt_orig; it++)
    if (TIG < 0 || int(it) == TIG) 
      ValidateEfastaRecord(tigs_orig[it]);


  // Load markups.

  vec< triple<int,int,int> > markups;
     
  if (USE_MARKUPS)
    markups_read(SCAFFOLDS_IN, & markups, tigs_orig, TIG, VERBOSE);



  // Go through the contigs.


  cout << Date() << ": finding edits" << endl;
  vec< vec<Edit> > edits_vec(nt_orig);
  size_t n_edits_tot = 0;

  for (size_t it = 0; it < nt_orig; it++) {
    if (TIG < 0 || int(it) == TIG) {

      const efasta & tig_orig = tigs_orig[it];

      vec<Edit> edits;
      edits_compile(tig_orig, markups, it, & edits);

      // Identify ambiguities that are tandem repeats having sufficiently large
      // period and sufficiently large ambiguity-free flanks.  
      // Form revised list of edits.
       
      const size_t n_edits = edits.size();
      for (size_t ie = 0; ie < n_edits; ie++) {

        Edit & edit = edits[ie];
        ValidateEfastaRecord(edit.efa);
           
        const bool valid = handle_true_tandem(tig_orig, &edit, MIN_PERIOD, USE_TANDEM);

        if (valid) {
          const int i0_flank = edit.i0 - FLANK;
          const int i1_flank = edit.i1 + FLANK;

          // make sure edits don't overlap each other

          if (i0_flank > 0 && 
              i1_flank < tig_orig.isize() &&
              (ie == 0           || i0_flank > edits[ie-1].i1) &&
              (ie == n_edits - 1 || i1_flank < edits[ie+1].i0)) {
               
            const String left  = tig_orig.substr(i0_flank, FLANK);
            const String right = tig_orig.substr(edit.i1, FLANK);
            const String all   = tig_orig.substr(i0_flank, i1_flank - i0_flank);
               
            if (VERBOSE) PRINT7(it, i0_flank, i1_flank, left, edit.efa, right, all);
               
            edits_vec[it].push(i0_flank, i1_flank, left + edit.efa + right);
            n_edits_tot++;
          }

        }

      }
       
    }
  }
  cout << Date() << ": found " << n_edits_tot << " edits" << endl;









  // Go through the scaffolds and edit them.  
  // We set the deviation of the new gaps to zero.


  cout << Date() << ": modifying scaffolds" << endl;
  vec<efasta> tigs_edit;
  size_t nt_tot = 0;

  vec<assembly_edit> a_edit2;
  vec<int> i2;

  const size_t ns = scaffolds.size();
  for (size_t is = 0; is < ns; is++) {

    superb & scaffold = scaffolds[is];
       
    for (int ist = scaffold.Ntigs() - 1; ist >= 0; ist--) {
      size_t it = scaffold.Tig(ist);

      const efasta & tig_orig = tigs_orig[it];

      const vec<Edit> & edits = edits_vec[it];
      const size_t n_edits = edits.size();
         
      const size_t nt_edit_loc = n_edits + 1;
      vec<efasta> tigs_edit_loc(nt_edit_loc);

      if (edits.empty()) {
        tigs_edit_loc[0] = tig_orig;
      }
      else {
        tigs_edit_loc[0] = tig_orig.substr(0, edits[0].i0);
        for (size_t l = 0; l < n_edits - 1; l++)
          tigs_edit_loc[l+1] = tig_orig.substr(edits[l].i1, 
                                               edits[l+1].i0 - edits[l].i1);
         
        tigs_edit_loc.back() = tig_orig.substr(edits.back().i1, 
                                               tig_orig.isize() - edits.back().i1);
      }



      vec<int> gap2(n_edits);

      for (size_t l = 0; l < n_edits; l++) {
        gap2[l] = edits[l].efa.Length1();
        vec<basevector> bvv_rep;
        edits[l].efa.ExpandTo(bvv_rep);
        a_edit2.push(assembly_edit::GAP_CLOSER, nt_tot + l, 
                     tigs_edit_loc[l].Length1(), nt_tot + l + 1, 0, bvv_rep);
        i2.push_back(nt_tot + l);
      }




      superb scaffold_edit;
      scaffold_edit.SetNtigs(nt_edit_loc);

      for (size_t itl = 0; itl < nt_edit_loc; itl++) {

        scaffold_edit.SetTig(itl, nt_tot + itl);
        scaffold_edit.SetLen(itl, tigs_edit_loc[itl].Length1());

        if (int(itl) < scaffold_edit.Ngaps()) {
          scaffold_edit.SetGap(itl, gap2[itl]);
          scaffold_edit.SetDev(itl, 0);
        }
      }



      // update scaffold and new contigs

      scaffold.ReplaceTigBySuper(ist, scaffold_edit);
         
      for (size_t it2 = 0; it2 < tigs_edit_loc.size(); it2++)
        tigs_edit.push_back(tigs_edit_loc[it2]);

      nt_tot = tigs_edit.size();

      i2.push_back(nt_tot - 1);
    }
  }







  // Write assembly and edits.

  if (WRITE) {
    cout << Date() << ": writing assembly and edits" << endl;
    Assembly A(scaffolds, tigs_edit);
    A.check_integrity();
    A.WriteAll(SCAFFOLDS_OUT);
    String outputFile = SCAFFOLDS_OUT + ".edits";
    BinaryWriter::writeFile(outputFile.c_str(), a_edit2);    
    UniqueSort(i2);
    Ofstream(out, SCAFFOLDS_OUT + ".tig_ids");
    for (int j = 0; j < i2.isize(); j++)
      out << i2[j] << "\n";
  }
}
