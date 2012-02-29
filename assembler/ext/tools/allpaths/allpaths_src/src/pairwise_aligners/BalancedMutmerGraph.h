/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2009) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

#ifndef BALANCED_MUTMERGRAPH
#define BALANCED_MUTMERGRAPH

#include "Basevector.h"
#include "CoreTools.h"
#include "math/Functions.h"
#include "pairwise_aligners/MaxMutmerFromMer.h"
#include "system/LockedData.h"

// ================================================================================
//
// A BMG (Balanced Mutmer Graph) is a vector of BMG_node's.  Each node has an
// array of mutmer_read_id's, stored as a vec; it might be possible to optimize for
// speed the way in which memory is allocated for the mutmer_read_ids.
//
// The mutmer_read_ids come in small and big versions,
// depending on the the value of a template parameter I (1 or 2).  If I = 1, then
// the mutmer_read_id occupies 8 bytes, and contains the following data:
// 
// * read id of child (31 bits)
// * reverse complement? (1 bit)
// * pos1 = position on parent read (10 bits)
// * pos2 = position on child read (10 bits)
// * length = length of overlap (10 bits)
// * unused bits (2 bits)
//
// If I = 2, you get the big version, which occupies 12 bytes:
//
// * read id of child (47 bits)
// * reverse complement? (1 bit)
// * pos1 = position on parent read (16 bits)
// * pos2 = position on child read (16 bits)
// * length = length of overlap (16 bits)
//
// The maximum dataset sizes are:
// I = 1: 2^31 reads, max read length 2^10
// I = 2: 2^47 reads, max read length 2^16
//
// The read id of the the parent is the index in the BMG of the
// BMG_node.
//
// If the reverse complement bit is set, then the region of overlap agrees as a
// reverse complement, e.g. as in:
//
//       ........ACGTT...    read 1      
//          .....AACGT.....  read 2
//
// This file should work for little endian architectures, but will need
// rewriting for the big endian case.
//
// ================================================================================

#ifndef Little_Endian
     #error BalancedMutmerGraph.h is designed for little endian architectures.
#endif





template<int I = 1> class mutmer_read_id 
{
private:
  uint32_t data_[1+I]; // i.e. 2*4 bytes if I=1, 3*4 bytes if I=2
  
  static const int MAX_READ_LEN_1 = 1024; // Maximum read length for I=1: 2^10

public:
  
  // Implicit assumptions (not checked, to save time):
  // I == 1 or 2
  // If I == 1, then read_id < 2^31, and pos1,pos2,len < 1024
  void Pack(const int64_t & read_id, 
            const uint16_t & pos1,
            const uint16_t & pos2, 
            const uint16_t & len,
            const Bool & rc) 
  {
      // Split the 31 least-significant bits off of read_id.
    // If I=1, these should be the ONLY bits of read_id.

    uint32_t read_id_0 = unsigned( read_id & (int64_t)0x7fffffff );

    // Byte-pack data_[0] with the id and RC flag.

    data_[0] = ( read_id_0 << 1 ) ^ rc;
       
    // Byte-pack data_[1] (and maybe data_[2]) with position information.
    // The bit-packing here mirrors that of the Unpack() function, below.

    if ( I == 1 ) {
      //ForceAssertLt( pos1, MAX_READ_LEN_1 );
      //ForceAssertLt( pos2, MAX_READ_LEN_1 );
      //ForceAssertLt( len,  MAX_READ_LEN_1 );
    
      data_[1] = pos1 | (pos2 << 10) | (len << 20);
    }
    else if ( I == 2 ) { // For I=2, get the next 16 bits from read_id.
      
      uint32_t read_id_31 = unsigned( read_id >> 31 ) & 0xffff;
      
      data_[1] = read_id_31 | (len << 16);
      data_[2] = pos1 | (pos2 << 16);
    }
  }

  // The bit-packing here mirrors that of the Pack() function, above.
  void Unpack(uint16_t & pos1, 
              uint16_t & pos2, 
              uint16_t & len) const 
  {
    if ( I == 1 ) {
      pos1 = data_[1] & 01777;
      pos2 = (data_[1] >> 10) & 01777;
      len = (data_[1] >> 20) & 01777;
    }
    else if ( I == 2 ) {
      pos1 = data_[2] & 0xffff;
      pos2 = (data_[2] >> 16) & 0xffff;
      len = (data_[1] >> 16) & 0xffff;
    }
  }


  int64_t UnpackReadId() const 
  {
    if ( I == 1 ) return data_[0] >> 1;
    else { // I=2
      int64_t read_id = data_[0] >> 1;
      read_id |= ( ( int64_t ( data_[1] & 0xffff ) ) << 31 );
      return read_id;
    }
  }
  
  inline Bool Rc() const 
  { return data_[0] & 1; }

  friend Bool operator==( const mutmer_read_id<I>& m1, const mutmer_read_id<I>& m2 )
  { return m1.data_[0] == m2.data_[0]; }

  friend Bool operator<(const mutmer_read_id<I>& m1, const mutmer_read_id<I>& m2)
  { return m1.data_[0] < m2.data_[0]; }
  
};




// Node in a Balanced Mutmer Graph
template<int I> 
class BMG_node : public vec< mutmer_read_id<I> > 
{
public:

  size_t SizeOf() const { return (*this).size() * sizeof(mutmer_read_id<I>); } 




  // ---- input from a stream
  friend std::istream& operator>>(std::istream & is, BMG_node & nodes) 
  {
    uint32_t size; // stored as 32 bits on disk
    is.read(reinterpret_cast<char*>(&size), sizeof(uint32_t));
    nodes.resize(size);

    //static int nn = 0;
    //if (nn++ < 10) cout << size << endl;

    for (size_t i = 0; i != size; i++)
      is.read(reinterpret_cast<char*>(&(nodes[i])), sizeof(mutmer_read_id<I>));

    return is; 
  }
  


  // ---- output to a stream
  friend std::ostream& operator<<(std::ostream & os, const BMG_node & nodes)
  { 
    uint32_t size = nodes.size(); // stored as 32 bits on disk
    os.write(reinterpret_cast<const char*>(&size), sizeof(uint32_t));
    
    //static int nn = 0;
    //if (nn++ < 10) cout << size << endl;

    for (size_t i = 0; i != size; i++)
      os.write(reinterpret_cast<const char*>(&(nodes[i])), sizeof(mutmer_read_id<I>));
    
    return os; 
  }
};












// Balanced Mutmer Graph
// Name is now a misnomer:  No balancing occurs.  The mutmer_read_id is always
// added to readID2's vector.  MergeKmer is now thread-safe.
template<int I> class BMG : public vec< BMG_node<I> > 
{
public:
  size_t SizeOf() const 
  {
    size_t bytes = 0;
    for (size_t i = 0; i < (*this).size(); i++) 
      bytes += (*this)[i].SizeOf();
    return bytes;
  }

  inline void push_back(const pair< int, mutmer_read_id<I> > & merge_node) 
  {   
    if (merge_node.first != -1)
      (*this)[merge_node.first].push_back(merge_node.second);
  }

  void FromFile(const String & fn)
  {
    ifstream is;
    is.open(fn.c_str());
    uint64_t n; // fixed at 64 bits on file
    is.read(reinterpret_cast<char*>(&n), sizeof(uint64_t));
    //cout << "FromFile n = " << n << endl;
    (*this).resize(n);
    for (size_t i = 0; i != n; i++)
      is >> (*this)[i];
    is.close();
  }

  void ToFile(const String & fn) const 
  {
    ofstream os;
    os.open(fn.c_str());
    uint64_t n = (*this).size();  // fixed at 64 bits on file
    //cout << "ToFile n = " << n << endl;
    os.write(reinterpret_cast<const char*>(&(n)), sizeof(uint64_t));
    for (size_t i = 0; i != n; i++)
      os << (*this)[i];
    os.close();
  }


  


  
  // MergeKmer: now a wrapper for MergeKmerConst
  void MergeKmer( int pos1, 
                  int pos2, 
                  int K, 
                  const int64_t& read_id1, 
                  const int64_t& read_id2, 
                  const basevector& rd1, 
                  const basevector& rd2,
                  Bool strict = False ) 
  {
    // Set rc if exactly one of pos1, pos2 are negative.  Then 
    // set pos_i to abs(pos_i) - 1, for i = 1, 2.
    
    Bool rc = (pos1 < 0) ^ (pos2 < 0);
    if ( pos1 < 0 ) pos1 = -pos1;
    if ( pos2 < 0 ) pos2 = -pos2;
    --pos1;
    --pos2;
    if ( rc )
    {
       pos1 = rd1.size() - K - pos1; // flip
    }
    

    //ForceAssert( pos1 >= 0 && pos2 >= 0 );
    //ForceAssert( pos1 < 65536 && pos2 < 65536 );
    
    uint16_t xpos1, xpos2, xk;

    // mutex protecting all nodes sharing lowest 10 bits of their ids
    Locker locker( mLocks[read_id2&0x3FF] );
    
    // Search the list associated with read_id2
    BMG_node<I>& node = (*this)[read_id2];
    typedef typename BMG_node<I>::iterator NodeItr;
    NodeItr end(node.end());
    for ( NodeItr itr(node.begin()); itr != end; ++itr )
    {
      if (itr->UnpackReadId() == read_id1 && rc == itr->Rc()) {
        itr->Unpack(xpos1, xpos2, xk);
	if (xpos1 - xpos2 == pos2 - pos1 &&
	    pos2          >= xpos1 &&
	    pos2 + K      <= xpos1 + xk)
	  return;
      }
    }

    int scratch = 0;
    if (rc) {
      MaxMutmerFromMerRev(pos2, pos1, K, scratch, rd2, rd1, strict);
    }
    else {
      MaxMutmerFromMer(pos2, pos1, K, scratch, rd2, rd1, strict);
    }

    mutmer_read_id<I> id;
    id.Pack(read_id1, pos2, pos1, K, rc);
    node.push_back(id);
  }

private:
  LockedData mLocks[1024];
};

#endif
