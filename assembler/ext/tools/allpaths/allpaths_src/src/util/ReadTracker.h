/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2008) by the   //
//   Broad Institute.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////


/* ReadTracker
 *
 * A quick-and-dirty class to manage tracking of read numbering
 * between modules which filter or combine reads.
 * 
 * Basic use by a program which geneates new reads file:
 *
 * ReadTracker rt;
 *
 * unsigned int source1 = rt.AddSource(source_filename1);
 * unsigned int source2 = rt.AddSource(source_filename2);
 * ...
 * 
 * for (<reads in new order>) {
 *   rt.AddRead(sourceX, sourceX_read_id);
 * }
 * 
 * rt.Dump(new_read_file_head);
 *
 * Programs which want to reference a previously generated ReadTracker file:
 *
 * ReadTracker rt;
 *
 * rt.Load(read_file_head);
 * String source = rt.GetReadSource(read_id);
 * uint64_t source_read_id = rt.GetReadSource(read_id);
 *
 * Bruce Walker
 * 16 Dec 09
 ******************************************************************************/

#define READTRACKER_EXT		".readtrack"
#define READTRACKER_SOURCE_BITS 8
#define READTRACKER_READID_BITS (64-READTRACKER_SOURCE_BITS)

#define RT_BITMASK(n) (((uint64_t)1<<(n))-1)

class ReadTracker {
public:
  vec<String> source_files;
  // read_index holds both source index and read_id in bitfields
  // expect this to change
  vec<uint64_t> read_index;

  unsigned int AddSource(String s);
  void AddReadSet( const String & source, const vec<Bool>& flags, const Bool keep_if_true = false );
  void AddRead(unsigned int source_index, uint64_t read_id);
  
  String GetSource(unsigned int source_index) const;
  unsigned int GetReadSourceIndex(unsigned int read) const;
  String GetReadSource(unsigned int read) const;
  uint64_t GetReadIndex(unsigned int read) const;
  void Dump(String filename) const;
  void Load(String filename);
  uint64_t size() const;
};
