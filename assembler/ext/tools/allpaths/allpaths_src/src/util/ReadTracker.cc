/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2009) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

// Interface documented (barely) in Readtracker.h

#include "VecUtilities.h"
#include "util/ReadTracker.h"

unsigned int 
ReadTracker::AddSource(String s)
{
  unsigned int n = source_files.size();
  source_files.push_back(s);
  return n;
}

void 
ReadTracker::AddRead(unsigned int source, uint64_t read_id)
{
  ForceAssertLt(source, 1u<<READTRACKER_SOURCE_BITS);
  ForceAssertLt(read_id, 1ul<<READTRACKER_READID_BITS);
  uint64_t x = read_id | ((uint64_t)source<<READTRACKER_READID_BITS);
  read_index.push_back(x);
}

// Add a whole set of reads at once.
// The flags indicate which reads to discard (or which to keep.)
void
ReadTracker::AddReadSet( const String & source, const vec<Bool> & flags,
			 const Bool keep_if_true )
{
  uint n = AddSource( source );
  for ( size_t i = 0; i < flags.size(); i++ )
    if ( flags[i] == keep_if_true )
      AddRead( n, i );
}


String
ReadTracker::GetSource(unsigned int source_index) const
{
  return source_files[source_index];
}


unsigned int
ReadTracker::GetReadSourceIndex(unsigned int read) const
{
  return (read_index[read] >> READTRACKER_READID_BITS)
    & RT_BITMASK(READTRACKER_SOURCE_BITS);
}

String
ReadTracker::GetReadSource(unsigned int read) const
{
  return GetSource(GetReadSourceIndex(read));
}


uint64_t 
ReadTracker::GetReadIndex(unsigned int read) const
{
  return read_index[read] & RT_BITMASK(READTRACKER_READID_BITS);
}

void
ReadTracker::Dump(String filename) const
{
  String filename_ext = filename + READTRACKER_EXT;
  ofstream f(filename_ext.c_str());
  for (size_t i = 0; i < source_files.size(); ++i)
    f << "S " << source_files[i] << "\n";
  for (size_t i = 0; i < read_index.size(); ++i)
    f << GetReadSourceIndex(i) << " " << GetReadIndex(i) << "\n";
  f.close();
}

void
ReadTracker::Load(String filename)
{
  source_files.clear();
  read_index.clear();
  String filename_ext = filename + READTRACKER_EXT;
  ifstream f(filename_ext.c_str());
  if (!f.good()) return;
  // read the header
  bool header = true;
  while(header){
    char c = f.peek();
    if ( c == 'S' ){
      string line;
      getline(f,line);
      istringstream iss(line);
      string token, name;
      iss>> token >> name;
      //cout << name << endl;
      AddSource(name);
    }
    else{
      header = false;
    }
  }
  // read the records
  uint32_t source=0;
  uint64_t read=0;
  long int i=0;
  const int bufSize = 10*1024*1024; // 10M reading buffer
  const int safeSize = 100; // additional 100 char safe zone 
  char* buf = new char[bufSize]; 
  // determine the start and end position of the records in the file
  long pos_start = f.tellg();
  f.seekg(0,ios::end);
  long pos_end = f.tellg();
  f.seekg(pos_start);
  while(1){
    long fstart = f.tellg();
    f.read(buf, bufSize);
    long pos = f.tellg();
    // get the read_length
    if (f.eof()) pos = pos_end;
    long read_length = pos - fstart;
    char* pLimit = buf + read_length - safeSize; // no converting beyond the safe zone
    char* p =  buf;
    while( p < pLimit) {
      source = strtoul(p,&p,10);
      read  = strtoul(p,&p,10);
      AddRead(source,read);	
      i++;
    } 
    // reload the string in the buffer not read
    long actual_read_length = p - buf;
    if ( f.eof()) {
      f.clear();
      f.seekg(fstart + actual_read_length);
      break;
    }else{
      f.seekg(fstart + actual_read_length);
    }
  }
  // the rest 
  while(f>>source)
  {
    f >> read;
    AddRead(source,read);	
    i++;
  }
  delete [] buf;
  //cout << "total reads " << i << endl;
  //cout << "last read " << source << " "<<read << endl;
  f.close();
}

uint64_t
ReadTracker::size() const
{
  return read_index.size();
}
