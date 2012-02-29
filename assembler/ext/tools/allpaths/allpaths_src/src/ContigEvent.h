// Copyright (c) 2004 Broad Institute/Massachusetts Institute of Technology

// NO LONGER SUPPORTED.

// This file defines classes contig_event and contig_event_list, which 
// facilitate the tracking of changes to contigs in MergeSupercontigs
// (and could be used elsewhere as well).  See also MergeSupercontigs.cc and
// TraceContig.cc.
//
// How to generate a contig_event_list (example):
//      Ofstream( out, "contig_event_file" );
//      {    merge_contig_event e(3, 4, 3);
//           e.BinaryWrite(out);    }
//      {    reverse_contig_event e(3);
//           e.BinaryWrite(out);    }
//      {    rename_contig_event e(3, 2);
//           e.BinaryWrite(out);    }

// Note.  There is a design defect in this code: we assume that in a derived class
// of an abstract class, the first item stored is the virtual table pointer.
// While under gcc this appears to always be the case, it may not be the case under
// other compilers.  The code does do a sanity test (SanityTest) -- if the 
// assumption made here is wrong, then an appropriate fatal error is issued.

#ifndef CONTIG_EVENT
#define CONTIG_EVENT

#include "system/Assert.h"
#include "String.h"
#include "system/System.h"
#include "system/Types.h"
#include "Vec.h"

class int_pointer_or_shorts {
     public:
     union 
     {    int* p;
          unsigned short x[2];    };
};

// ===========================================================================

class contig_event {

     public:

     static vec<String> class_names;
     static vec<void*> virtual_table_pointers;

     // size of an event in bytes, or in machine words

     virtual int SizeOf( ) const = 0;
     int SizeOfInWords( ) const
     {    int n = sizeof(int*);
          return (SizeOf( ) + n - 1) / n;    }

     virtual char* Loc( ) const = 0;

     virtual vec<int> InputContigs( ) const = 0;
     virtual int OutputContig( ) const = 0;

     virtual Bool Trivial( ) const = 0;
     virtual Bool Marker( ) const = 0;
     virtual Bool ContigOperator( ) const = 0;
     virtual Bool LinkCreationOperator( ) const = 0;
     virtual Bool ReverseLinkOperator( ) const = 0;

     virtual String ContigEventClassName( ) const = 0;
     int ContigEventClassId( ) const;

     // human-readable output, or binary output

     virtual void Print( ostream& out ) const = 0;

     void BinaryWrite( ostream& out ) const
     {    ForceAssert( sizeof(int_pointer_or_shorts) == sizeof(void*) );
          if ( Trivial( ) ) return;
          int_pointer_or_shorts p;
          unsigned short n = SizeOf( );
          p.x[0] = n;
          p.x[1] = ContigEventClassId( );
          BinWrite( out, p );
          out.write( Loc( ), n );    }

     Bool IsInput( int id ) const
     {    return Member( InputContigs( ), id );    }

     Bool IsOutputContig( int id ) const
     {    return id == OutputContig( );    }

public:
  virtual ~contig_event(){}
};

void InitializeContigEventStaticMembers( );

// ================================================================================

#define contig_event_class_prelude(EVENT_CLASS_NAME)                      \
                                                                          \
     class EVENT_CLASS_NAME : public contig_event {                       \
                                                                          \
     public:                                                              \
                                                                          \
     EVENT_CLASS_NAME( ) { }                                              \
     virtual ~EVENT_CLASS_NAME() {}                                       \
                                                                          \
     int SizeOf( ) const { return sizeof(*this) - sizeof(void*); }        \
     char* Loc( ) const { return (char*) this + sizeof(void*); }          \
                                                                          \
     String ContigEventClassName( ) const                                 \
     {    return #EVENT_CLASS_NAME;    }

// ================================================================================

contig_event_class_prelude(reverse_contig_event)

     reverse_contig_event( int id ) : id_(id) { }

     Bool ContigOperator( ) const { return True; }
     Bool LinkCreationOperator( ) const { return False; }
     Bool ReverseLinkOperator( ) const { return False; }

     vec<int> InputContigs( ) const
     {    vec<int> v(1);
          v[0] = id_;
          return v;    }

     int OutputContig( ) const { return id_; }

     Bool Trivial( ) const { return False; }
     Bool Marker( ) const { return False; }

     void Print( ostream& out ) const
     {    out << "reverse contig " << id_ << "\n";    }

     private: int id_;

};

// ================================================================================

contig_event_class_prelude(rename_contig_event)

     rename_contig_event( int id_old, int id_new ) 
          : id_old_(id_old), id_new_(id_new) { }

     Bool ContigOperator( ) const { return True; }
     Bool LinkCreationOperator( ) const { return False; }
     Bool ReverseLinkOperator( ) const { return False; }

     vec<int> InputContigs( ) const
     {    vec<int> v(1);
          v[0] = id_old_;
          return v;    }

     int OutputContig( ) const { return id_new_; }

     Bool Trivial( ) const { return id_old_ == id_new_; }
     Bool Marker( ) const { return False; }

     void Print( ostream& out ) const
     {    out << "rename contig " << id_old_ << " to " << id_new_ << "\n";    }

     // TODO: potentially dangerous truncation of index
     private: int id_old_, id_new_;

};


// ================================================================================

contig_event_class_prelude(merge_contig_event)

     merge_contig_event( int id_old1, int id_old2, int id_new ) 
          : id_old1_(id_old1), id_old2_(id_old2), id_new_(id_new) { }

     Bool ContigOperator( ) const { return True; }
     Bool LinkCreationOperator( ) const { return False; }
     Bool ReverseLinkOperator( ) const { return False; }

     vec<int> InputContigs( ) const
     {    vec<int> v(2);
          v[0] = id_old1_;
          v[1] = id_old2_;
          return v;    }

     int OutputContig( ) const { return id_new_; }

     Bool Trivial( ) const { return False; }
     Bool Marker( ) const { return False; }

     void Print( ostream& out ) const
     {    out << "merge contigs " << id_old1_ << " and " << id_old2_ << " to yield " 
               << id_new_ << "\n";    }

     private: int id_old1_, id_old2_, id_new_;

};

// ================================================================================

contig_event_class_prelude(merge_via_chimera_contig_event)

     merge_via_chimera_contig_event( int id_old1, int id_old2, int id_new ) 
          : id_old1_(id_old1), id_old2_(id_old2), id_new_(id_new) { }

     Bool ContigOperator( ) const { return True; }
     Bool LinkCreationOperator( ) const { return False; }
     Bool ReverseLinkOperator( ) const { return False; }

     vec<int> InputContigs( ) const
     {    vec<int> v(2);
          v[0] = id_old1_;
          v[1] = id_old2_;
          return v;    }

     int OutputContig( ) const { return id_new_; }

     Bool Trivial( ) const { return False; }
     Bool Marker( ) const { return False; }

     void Print( ostream& out ) const
     {    out << "merge contigs " << id_old1_ << " and " << id_old2_ 
               << " via a chimeric read to yield " << id_new_ << "\n";    }

     private: int id_old1_, id_old2_, id_new_;

};

// ================================================================================

contig_event_class_prelude(delete_contig_event)

     delete_contig_event( int id ) : id_(id) { }

     Bool ContigOperator( ) const { return True; }
     Bool LinkCreationOperator( ) const { return False; }
     Bool ReverseLinkOperator( ) const { return False; }

     vec<int> InputContigs( ) const
     {    vec<int> v(1);
          v[0] = id_;
          return v;    }

     int OutputContig( ) const { return -1; }

     Bool Trivial( ) const { return False; }
     Bool Marker( ) const { return False; }

     void Print( ostream& out ) const
     {    out << "delete contig " << id_ << "\n";    }

     private: int id_;

};

// ================================================================================

contig_event_class_prelude(link_contig_event)

     link_contig_event( int id1, int id2, int category ) 
          : id1_(id1), id2_(id2), category_(category) { }

     Bool ContigOperator( ) const { return False; }
     Bool LinkCreationOperator( ) const { return True; }
     Bool ReverseLinkOperator( ) const { return False; }

     vec<int> InputContigs( ) const
     {    vec<int> v(2);
          v[0] = id1_;
          v[1] = id2_;
          return v;    }

     int OutputContig( ) const { return -1; }

     Bool Trivial( ) const { return False; }
     Bool Marker( ) const { return False; }

     void Print( ostream& out ) const
     {    out << "link contigs " << id1_ << " and " << id2_ 
               << " (category = " << category_ << ")\n";    }

     private: int id1_, id2_, category_;

};

// ================================================================================

contig_event_class_prelude(reverse_link_contig_event)

     reverse_link_contig_event( int id1, int id2 ) : id1_(id1), id2_(id2) { }

     Bool ContigOperator( ) const { return False; }
     Bool LinkCreationOperator( ) const { return False; }
     Bool ReverseLinkOperator( ) const { return True; }

     vec<int> InputContigs( ) const
     {    vec<int> v(2);
          v[0] = id1_;
          v[1] = id2_;
          return v;    }

     int OutputContig( ) const { return -1; }

     Bool Trivial( ) const { return False; }
     Bool Marker( ) const { return False; }

     void Print( ostream& out ) const
     {    out << "reverse link between contigs " << id1_ << " and " << id2_ 
               << "\n";    }

     private: int id1_, id2_;

};

// ================================================================================

contig_event_class_prelude(mark_mergecontigs_start_contig_event)

     Bool ContigOperator( ) const { return False; }
     Bool LinkCreationOperator( ) const { return False; }
     Bool ReverseLinkOperator( ) const { return False; }
     vec<int> InputContigs( ) const { return vec<int>( ); }
     int OutputContig( ) const { return -1; }
     Bool Trivial( ) const { return False; }
     Bool Marker( ) const { return True; }
     void Print( ostream& out ) const { out << "start MergeContigsWithGaps\n"; }

};

// ================================================================================

contig_event_class_prelude(mark_mergesupercontigs_start_contig_event)

     Bool ContigOperator( ) const { return False; }
     Bool LinkCreationOperator( ) const { return False; }
     Bool ReverseLinkOperator( ) const { return False; }
     vec<int> InputContigs( ) const { return vec<int>( ); }
     int OutputContig( ) const { return -1; }
     Bool Trivial( ) const { return False; }
     Bool Marker( ) const { return True; }
     void Print( ostream& out ) const { out << "start MergeSupercontigs\n"; }

};

// ================================================================================

class contig_event_list {

     public:

     contig_event& operator[ ]( int i ) const
     {    return *((contig_event*) data_ + event_locs_[i]);    }

     void PrintAll( ostream& out )
     {    for ( unsigned int i = 0; i < event_locs_.size( ); i++ )
               (*this)[i].Print(out);    }

     contig_event_list( ) { }

     void SanityTest( )
     {    reverse_contig_event e(10234);
          int e_contig = *( (int*) ( ((int**) &e) + 1 ) );
          if ( e_contig != 10234 )
               FatalErr( "In internal error has been detected.  This error is\n"
                    << "probably due to a design defect described in "
                    << "ContigEvent.h." );    }

     // Construct from a file, created by a succession of BinaryWrite's.

     contig_event_list( const String& filename )
     {    SanityTest( );
          InitializeContigEventStaticMembers( );
          longlong n = FileSize(filename);
          int k = n / sizeof(int*);
          data_ = new int*[k];
          Ifstream( in, filename );
          in.read( (char*) data_, n );

          // Compute number of events.

          int event_count = 0, ptr = 0, nwords;
          int_pointer_or_shorts p;
          while( ptr < k )
          {    p.p = data_[ptr];
               nwords = 1 + (p.x[0] + sizeof(int*) - 1) / sizeof(int*);
               ptr += nwords;
               ++event_count;    }

          // Define event_locs, fill in table pointers.

          event_locs_.resize(event_count);
          ptr = 0;
          event_count = 0;
          unsigned short id;
          while( ptr < k )
          {    event_locs_[event_count] = ptr;
               p.p = data_[ptr];
               nwords = 1 + (p.x[0] + sizeof(int*) - 1) / sizeof(int*);
               id = p.x[1];
               data_[ptr] = (int*) contig_event::virtual_table_pointers[id];
               ptr += nwords;
               ++event_count;    }    }

     ~contig_event_list( ) { delete[ ] data_; }

     void Dump( );
     void TraceBack( int e );
     void TraceBack( int e1, int e2 );

     private:

     int** data_;
     vec<int> event_locs_;

};

#endif
