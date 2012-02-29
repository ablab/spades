///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////


#include "layout/ContigActualloc.h"

ostream& operator<<( ostream &o, const contig_actualloc &a ) {
  o << a.Weight() << " " 
    << a.Start()  << " "
    << BoolToInt( a.rc_ );
  return o;
}

istream& operator>>( istream &in, contig_actualloc &a ) {
  in >> a.weight_
     >> a.start_;
  int dummy;
  in >> dummy;
  a.rc_ = ( Bool ) dummy;

  return in;
}


bool operator<( const contig_actualloc &a1, const contig_actualloc &a2 ) {
  if ( a1.Weight() != a2.Weight() ) return a1.Weight() < a2.Weight();
  else                              return a1.Start()  < a2.Start();
}



ostream& operator<<( ostream &o, const arachne_contig &a ) {
  o << a.id_     << " "
    << a.length_ << " "
    << a.sc_id_  << " "
    << a.sc_pos_ << endl;

  o << a.actuallocs_.size() << endl;
  ctg_actloc_itr curr = a.actuallocs_.begin();
  ctg_actloc_itr end  = a.actuallocs_.end();

  for ( ;
	curr != end;
	++curr )
    o << *curr << endl;

  return o;
}

istream& operator>>( istream &in, arachne_contig &a ) {
  in >> a.id_
     >> a.length_
     >> a.sc_id_
     >> a.sc_pos_;
  
  int len;
  in >> len;
  Assert( len >= 0 );

  for ( int i = 0;
	i < len;
	++i ) {
    contig_actualloc cactl;
    in >> cactl;
    a.actuallocs_.insert( cactl );
  }

  return in;
}

void contig_actualloc::Print( ostream &o ) const {
  o << "<weight: " << weight_ 
    << "; start,dir: " << start_ 
    << "." << (int)rc_ 
    << ">";
}


void arachne_contig::Print( ostream &o ) const {
  o << "ID: " << id_ << "; SC = " << sc_id_ << ":" << sc_pos_ 
    << "; l = " << length_
    << "; actlocs: { ";

  ctg_actloc_itr curractloc_itr = actuallocs_.begin();
  ctg_actloc_itr end_actloc_itr = actuallocs_.end();

  for ( ;
        curractloc_itr != end_actloc_itr;
        ++curractloc_itr ) {
    curractloc_itr->Print( o );
    o << ",";
  }
  o << "}";
}


