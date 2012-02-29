/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2008) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

#include "CoreTools.h"
#include "graph/Digraph.h"
#include "graph/DigraphTemplate.h"
#include "paths/UnipathNhoodCommon.h"
#include "paths/simulation/Placement.h"


/**
   Function: PrintNhood

   Prints the <neighborhood> to stdout, for debugging purposes.

   Parameters:

      v - the <neighborhood seed>
      processed - the <neighborhood unipaths>
      ulen - the lengths of unipaths
      USE_TRUTH - whether we have the reference and have the placements of unipaths
         on the reference.
      pLocs - if USE_TRUTH, then for each unipath, the placements of that unipath
        on the reference genome.
*/
void PrintNhood( int v, const vec<ustart>& processed, const vec<nbases_t>& ulen,
     Bool USE_TRUTH, const VecPlacementVec* pLocs )
{    if ( USE_TRUTH ) { ForceAssert( pLocs != 0 ); }
     cout << "\nneighborhood of " << v << ":\n";
     vec< vec<String> > rows;
     bool seed_is_copy_one = ( USE_TRUTH && (*pLocs)[v].size() == 1 );
     for ( int j = 0; j < processed.isize( ); j++ )
     {    int start = processed[j].Start( );
          int w = processed[j].Uid( );
	  vec<String> row;
          ostrstream out0, out1, out2;
          out0 << j+1 << ends;
          out1 << "[" << start << "," << start + ulen[w] - 1 << "]"
               << " +/- " << processed[j].MeanDev( ) << ends; 
          out2 << w << " (" << ulen[w] << " kmers)" << ends;
          row.push_back( out0.str( ), out1.str( ), out2.str( ) );
          if (USE_TRUTH)
          {    ostrstream out3;
               for ( PlacementVec::size_type u = 0; u < (*pLocs)[w].size( ); u++ )
               {    if ( u > 0 ) out3 << ", ";
                    const placement& thisPlacement = (*pLocs)[w][u];
                    out3 << thisPlacement;
                    if ( w != v && seed_is_copy_one ) {
                      const placement& seedPlacement = (*pLocs)[v].front();
                      if ( thisPlacement.GenomeId() == seedPlacement.GenomeId() &&
                           thisPlacement.Rc() == seedPlacement.Rc() ) {
                        int actualStart = 
                          ( thisPlacement.Rc() ? 
                            seedPlacement.Pos() - thisPlacement.Pos() :
                            thisPlacement.pos() - seedPlacement.pos() );
                        int diff = start - actualStart;
                        float devsOff = (float)diff / (float)processed[j].MeanDev();
                        /*
                        out3 << " [diff=" << diff
                             << ", devs=" << setprecision(2) << devsOff << "]";
                        */
                      }
                    }
               }
               out3 << ends;
               row.push_back( out3.str( ) );
          }
          rows.push_back(row);    }
     PrintTabular( cout, rows, 1 );    }



