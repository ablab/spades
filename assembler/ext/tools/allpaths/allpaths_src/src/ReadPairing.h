///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////


#ifndef READTYPE
#define READTYPE

#include "String.h"
#include "system/Types.h"
#include "Vec.h"


enum read_type { insert_4k,
                 insert_40k,
                 gFgR_transposon,
                 g7g8_transposon,
                 G3_subclone,
                 preassembled_contig_piece,
                 other,
                 shattered_pcr_product,
		 bac_clone,
		 long_bac_clone,
		 plasmid_1000,
		 plasmid_2000,
		 plasmid_3000,
		 plasmid_4000,
		 plasmid_5000,
		 plasmid_6000,
		 plasmid_7000,
		 plasmid_8000,
		 plasmid_9000,
		 plasmid_10000,
		 unspecified_read_type
};

/*
   Class: read_pairing

   Represents a pairing of two reads, indicating that they come from two ends of
   the same DNA fragment whose length is roughly known (and so the distance
   between these two reads in the genome is roughly known also).

   The reads of the pair point _towards_ each  other:
   >    --------->.............................<---------
   >      read1     unread part of fragment       read2
*/
class read_pairing {

     public:

  // Fields: read identifiers
  //
  // Identifies the reads of the pair.  (Which read is first and which is second
  // is arbitrary).
  //
  //    id1 - <read identifier> of the first read
  //    id2 - <read identifier> of the second read
  // TODO: potentially dangerous truncation of index by id1 and id2
  int id1, id2;   // identifiers of the two reads

  // Fields: read separation
  //
  // Mean and standard deviation of the separation between the reads.  We don't
  // know the actual separation, unless we know the reference genome, so
  // normally we'll just know the mean and standard deviation of the fragment
  // sizes in a library; while here we're representing these values individually
  // for every read pair, normally the read separations all come from the same
  // distribution and will have the same values for all read pairings.
  // 
  //    sep - mean value of the separation between the reads
  //    sd - standard deviation of the separation between the reads
  int sep;        
  int sd;

  // Field: t
  //
  // <Read type> of the two reads.
  read_type t;

  // Field: weight
  //
  // Weight of the pairing (what is this -- how sure we are that the pairing is correct?)
  int weight;     

     read_pairing( ) { id1 = -1; }
     read_pairing( int _id1, int _id2, int _sep, int _sd ):
       id1(_id1), id2(_id2), sep(_sep), sd(_sd), t(unspecified_read_type), weight(0) { }
     read_pairing( int _id1, int _id2, int _sep, int _sd, read_type _t ):
       id1(_id1), id2(_id2), sep(_sep), sd(_sd), t(_t), weight(0) { }

     void ChangeId( int old_id, int new_id );

     bool InvolvesId( int id ) const;

     int Partner( int id ) const;

     void Kill( );

     bool Dead( ) const { return id1 == -1; }
     bool Alive( ) const { return id1 != -1; }

     // swaps the first read of the pair with the second.
     void Swap();

     bool operator==( const read_pairing &other ) const
     {
         return ( id1 == other.id1 &&
                  id2 == other.id2 &&
                  sep == other.sep &&
                  sd  == other.sd  &&
                  t   == other.t   &&
                  weight == other.weight );
     }

     friend Bool operator<( const read_pairing& p1, const read_pairing& p2 )
     {    return p1.id1 < p2.id1;    }

     friend ostream& operator<<( ostream& o, const read_pairing& p );
     friend istream& operator>>( istream& i, read_pairing& p );
     friend istream& operator>>( istream& i, vec<read_pairing>& v );

};

void ReadBinaryPairs ( istream &in, vec<read_pairing> &vecPairs );

// WritePairs: write files reads.pairto, reads.pairtob, reads.pairto_index.

void WritePairs( const String& dir, const vec<read_pairing>& pairs, size_t nreads,
		 Bool extras_only = False, const String& reads_name = "reads" );

void WritePairs( const vec<read_pairing>& pairs, size_t nreads, const String& file_head,
		 Bool extras_only = False );

// ReadsToPairs: given reads "ids", return the corresponding read pairs (leaving
// a dead pair in case the read is unpaired).

void ReadsToPairs( const String& dir, const vec<int>& ids,
		   vec<read_pairing>& pairs, const String& reads_name = "reads");

// If there is a binary file (strTextFilename+"b") newer than text file
// (strTextFilename), will read in binary file.  Otherwise, reads in text file.
void ReadPairsFile( const String &strTextFilename, vec<read_pairing> &vecPairs );


// If there is a binary file (strTextFilename+"b") newer than text file
// (strTextFilename), will read number of pairs from binary file.  Otherwise,
// read the number of pairs from the text file.
int PairsFileSize( const String &strTextFilename );


// Logical Type: pair_id_t
// Logical type for a read pair id -- the index of a read pair in vector of read pairs.
typedef int pair_id_t;

#ifdef __DECCXX_VER
#pragma do_not_instantiate istream& operator>>(istream&, vec<read_pairing>&)
#endif

#endif

// Synonyms: Various synonyms
//   read pairing - See <read_pairing>.
