///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// This code allows some tracking of news and deletes.

// We specifically do not include MemTracker.h here.

#include <stdlib.h>
#include <unistd.h>
#include <new>

#ifdef TRACK_MEMORY

#undef new

// As operator new() is called, we build up a hash table of linked
// lists of "new records", which store the location of the memory
// allocated, where the allocation was requested, and a pointer to the
// next "new record" in the list.

// When operator delete() is called, we search the hash table for a
// matching allocation.  If it doesn't exist, then we assume the
// memory was allocated using the built-in operator new(), and do
// nothing.  If it exists and it's already been flagged as deleted, we
// print a warning message and perform a traceback.  If it exists and
// it hasn't been flagged as deleted, we flag it.

// At the end of execution, we print out information on all the memory
// blocks that were not freed.


#define USE_BOOKENDS 1

// If USE_BOOKENDS is defined, each block handed out by new has 16
// extra bytes allocated, 8 before the block the user is aware of, and
// 8 after.  Each of these 8 bytes is set to a particular value so
// that we can check for buffer overruns when the block is freed.


// A struct to track allocations via operator new.

typedef
struct tag_NEWRECORD
{
  struct tag_NEWRECORD *mp_next;
  void                 *mp_data;
  const char           *msz_file;
  std::size_t           m_size;
  int                   m_line;
  int                   mb_deleted;
} NewRecord;


// A block of these structs, to reduce the number of mallocs we do for
// our tracking.  Each block points to the previously allocated block.
// The records in each block can belong to different hash buckets.

const int numRecordsPerBlock = 2000;

typedef 
struct tag_NEWRECORDBLOCK
{
  struct tag_NEWRECORDBLOCK *mp_next;
  NewRecord                  m_records[ numRecordsPerBlock ];
  int                        m_nextFreeRecordIndex;
} NewRecordBlock;


// Allocation records are hashed by the location of the allocated
// memory modulo some reasonably large prime number.
const unsigned long bigPrime = 999983l;

// The hash table of linked lists of allocation records (one for
// object allocations and one for array allocations).
static NewRecord** sp_newRecordHashTables[2] = { 0, 0 };

enum {
 objectsTable = 0,
 arraysTable  = 1
};

static char* sz_unknownFile = "unknown";

static const char* sz_bufferBookend = "deaddead";
static const int   bufferBookendSize = 8;

// The linked list of allocation record blocks.
static NewRecordBlock* sp_newRecordBlockListHead = 0;

// Are we tracking memory?  We turn this off when we do tracebacks.
static bool s_trackingOn = true;

// Count it up!
static long memory_used = 0;

#include "system/Assert.h"


// A function to find a record in a hash table.  Return 0 if not found.
NewRecord* find_record( int hashTableIndex, void *pointer )
{
  NewRecord **hashTable = sp_newRecordHashTables[ hashTableIndex ];

  // If there is no hash table, do nothing.
  if ( hashTable == 0 )
    return 0;

  // Find the head of the list in the appropriate bucket.
  NewRecord *p_thisRecord 
    = hashTable[ reinterpret_cast<unsigned long>(pointer) % bigPrime ];

  // Walk through the bucket looking for the pointer.
  while ( p_thisRecord != 0 &&
          p_thisRecord->mp_data != pointer )
    p_thisRecord = p_thisRecord->mp_next;

  return p_thisRecord;
}


// A function to add a record to the list.
void add_record( int hashTableIndex, 
                 void *pointer, std::size_t size, const char* file, const int line )
{
  NewRecord ** &hashTable = sp_newRecordHashTables[ hashTableIndex ];

  // Allocate and initialize the hash table if it doesn't exist.
  if ( hashTable == 0 )
  {
    hashTable = (NewRecord**) ::malloc( bigPrime * sizeof(NewRecord*) );

    if ( hashTable == 0 )
    {
      fprintf( stderr, "\nUnable to allocate memory to track allocation from %s:%d.\n",
               file, line );
      // We turn off tracking so we don't put ourselves in an infinite
      // loop if the traceback has a bad new/delete.
      s_trackingOn = false;
      TracebackThisProcess();
    }

    for ( unsigned long index = 0; index < bigPrime; ++index )
      hashTable[index] = 0;
  }

  // Recycle the record if it already exists.
  NewRecord *p_thisRecord = find_record( hashTableIndex, pointer );

  if ( p_thisRecord == 0 )
  {
    // If we don't have any blocks or the current block is full,
    // allocate and initialize a new one.
    if ( sp_newRecordBlockListHead == 0 ||
         sp_newRecordBlockListHead->m_nextFreeRecordIndex == numRecordsPerBlock )
    {
      NewRecordBlock *p_newBlock = (NewRecordBlock*) ::malloc( sizeof(NewRecordBlock) );
      
      if ( p_newBlock == 0 )
      {
        fprintf( stderr, "\nUnable to allocate memory to track allocation from %s:%d.\n",
                 file, line );
        // We turn off tracking so we don't put ourselves in an infinite
        // loop if the traceback has a bad new/delete.
        s_trackingOn = false;
        TracebackThisProcess();
      }
      
      p_newBlock->m_nextFreeRecordIndex = 0;
      p_newBlock->mp_next = sp_newRecordBlockListHead;
      sp_newRecordBlockListHead = p_newBlock;
    }
    
    // Get a pointer to the next free record in the current block.
    p_thisRecord =
      sp_newRecordBlockListHead->m_records + 
      sp_newRecordBlockListHead->m_nextFreeRecordIndex; 
  
    ++sp_newRecordBlockListHead->m_nextFreeRecordIndex;
  
    // Figure out which bucket this record should go in.
    NewRecord **pp_newRecordListHead 
      = hashTable + (reinterpret_cast<unsigned long>(pointer) % bigPrime);
  
    // Stick the record in its bucket.
    p_thisRecord->mp_next  = *pp_newRecordListHead;
    *pp_newRecordListHead   = p_thisRecord;

    p_thisRecord->mp_data  = pointer;
  }

  // Add to the count
  memory_used += size;

  // Set the info for the record.
  p_thisRecord->m_size   = size;
  p_thisRecord->msz_file = file;
  p_thisRecord->m_line   = line;
  p_thisRecord->mb_deleted = 0;
}


// A function to mark a record as deleted (it is kept for reference).
// Returns a pointer to the record found, or null if no record was found.
NewRecord* delete_record( int hashTableIndex, void *pointer )
{
  // Find the head of the list in the appropriate bucket.
  NewRecord *p_thisRecord = find_record( hashTableIndex, pointer );

  // If we don't find it, traceback.
  if ( p_thisRecord == 0 )
  {
    fprintf( stderr, "\nUnable to find allocation record for pointer %p!\n",
             pointer );
    // We turn off tracking so we don't put ourselves in an infinite
    // loop if the traceback has a bad new/delete.
    s_trackingOn = false;
    TracebackThisProcess();
  }

  // If it's been deleted already, traceback.
  if ( p_thisRecord->mb_deleted )
  {
    fprintf( stderr, "\nDouble deletion of block of size %ld allocated at %s:%d!\n",
             p_thisRecord->m_size, 
             p_thisRecord->msz_file,
             p_thisRecord->m_line );
    // We turn off tracking so we don't put ourselves in an infinite
    // loop if the traceback has a bad new/delete.
    s_trackingOn = false;
    TracebackThisProcess();
  }

  // Flag it as deleted.
  p_thisRecord->mb_deleted = 1;

  // Subtract it from the count
  memory_used -= p_thisRecord->m_size;

  return p_thisRecord;
}


void * 
alloc_memory( std::size_t size, const char *file, const int line )
{
  if ( USE_BOOKENDS )
    size += 2 * bufferBookendSize;

  char *pointer = (char*) ::malloc( size );
  if ( pointer == 0 )
  {
    printf( "\nCall to new at %s:%d failed.\n", file, line );
    s_trackingOn = false;
    TracebackThisProcess();
  }

  if ( USE_BOOKENDS )
  {
    // Set the bookends to trap overruns.
    memcpy( pointer, sz_bufferBookend, bufferBookendSize );
    memcpy( pointer + size - bufferBookendSize, sz_bufferBookend, bufferBookendSize );
    
    pointer += bufferBookendSize;
  }

  return pointer;
}

// return true if bookends were tampered with
bool
check_bookends( NewRecord *p_record )
{
  bool tamperedWith = false;

  // Check the bookends for overruns.
  char *p_headBookend = ((char*) p_record->mp_data) - bufferBookendSize;
  char *p_tailBookend = ((char*) p_record->mp_data) + p_record->m_size;
    
  if ( memcmp( p_headBookend, sz_bufferBookend, bufferBookendSize ) ) {
    printf( "\nThe block allocated at %s:%d has had\n", 
	    p_record->msz_file, p_record->m_line );
    printf( "the %d bytes before it tampered with.\n", bufferBookendSize );
    tamperedWith = true;
  }
  
  if ( memcmp( p_tailBookend, sz_bufferBookend, bufferBookendSize ) ) {
    printf( "\nThe block allocated at %s:%d has had\n", 
	    p_record->msz_file, p_record->m_line );
    printf( "the %d bytes after it tampered with.\n", bufferBookendSize );
    tamperedWith = true;
  }

  return tamperedWith;
}

void
free_memory( void *pointer, NewRecord *p_record )
{
  if ( USE_BOOKENDS && p_record != 0 )
  {
    if ( check_bookends( p_record ) ) {
      s_trackingOn = false;
      TracebackThisProcess();
    }
    ::free( ((char*) pointer) - bufferBookendSize );
  }
  else
  {
    ::free( pointer );
  }
}

// User-defined braindead versions of the new and delete operators.

void * 
operator new ( std::size_t size, const char *file, const int line ) throw()
{
  void *pointer = alloc_memory( size, file, line );
  if ( s_trackingOn )
    add_record( objectsTable, pointer, size, file, line );
  return pointer;
}
  
void * 
operator new[] ( std::size_t size, const char *file, const int line ) throw()
{
  void *pointer = alloc_memory( size, file, line );
  if ( s_trackingOn )
    add_record( arraysTable, pointer, size, file, line );
  return pointer;
}

void * 
operator new ( std::size_t size ) throw(std::bad_alloc)
{
  void *pointer = ::operator new ( size, sz_unknownFile, 0 );
  if ( pointer == 0 )
    throw std::bad_alloc();
  return pointer;
}
  
void * 
operator new[] ( std::size_t size ) throw(std::bad_alloc)
{
  void *pointer = ::operator new[] ( size, sz_unknownFile, 0 );
  if ( pointer == 0 )
    throw std::bad_alloc();
  return pointer;
}

void * 
operator new ( std::size_t size, const std::nothrow_t& ) throw()
{
  return ::operator new ( size, sz_unknownFile, 0 );
}
  
void * 
operator new[] ( std::size_t size, const std::nothrow_t& ) throw()
{
  return ::operator new[] ( size, sz_unknownFile, 0 );
}

void
operator delete ( void * __p ) throw()
{
  if ( s_trackingOn && __p != 0 )
  {
    NewRecord *p_record = find_record( arraysTable, __p );
    if ( p_record != 0 &&
         p_record->mb_deleted == false )
    {
      printf( "\nAttempted to delete (as an object) a pointer returned by array new.\n" );
      s_trackingOn = false;
      TracebackThisProcess();
    }
    p_record = delete_record( objectsTable, __p );
    free_memory( __p, p_record );
  }
}

void
operator delete[] ( void * __p ) throw()
{
  if ( s_trackingOn && __p != 0 )
  {
    NewRecord *p_record = find_record( objectsTable, __p );
    if ( p_record != 0 &&
         p_record->mb_deleted == false )
    {
      printf( "\nAttempted to delete[] a pointer returned by object new.\n" );
      s_trackingOn = false;
      TracebackThisProcess();
    }
    p_record = delete_record( arraysTable, __p );
    free_memory( __p, p_record );
  }
}

void
operator delete ( void * __p, const std::nothrow_t& ) throw()
{
  ::operator delete ( __p );
}

void
operator delete[]( void * __p, const std::nothrow_t& ) throw()
{
  ::operator delete[] ( __p );
}

// Check the bookends of all of the memory allocations we're tracking.
// If any have been corrupted, print them out, and then a list of all
// currently tracked memory allocations, which may help figure out
// where the clobbering is happening.

void
check_malloc_bookends()
{
  if ( ! USE_BOOKENDS ) {
    fprintf( stderr, "\nMalloc bookends are not in use, so they cannot be checked.\n" );
    return;
  }

  bool bookendTrampled = false;
  for ( int listIndex = 0; listIndex < 2; ++listIndex )
  {
    NewRecord** hashTable = sp_newRecordHashTables[listIndex];

    if ( hashTable != 0 )
    {
      for ( unsigned int index = 0; index < bigPrime; ++index )
      {
        NewRecord *p_thisRecord = hashTable[index];
        
        while ( p_thisRecord != 0 ) {
          if ( ! p_thisRecord->mb_deleted &&
               p_thisRecord->msz_file != sz_unknownFile )
	    if ( check_bookends( p_thisRecord ) ) 
	      bookendTrampled = true;
          p_thisRecord = p_thisRecord->mp_next;
	}
      }
    }
  }

  if ( bookendTrampled ) {
    fprintf( stderr, "\nList of all current memory allocations:\n" );
    for ( int listIndex = 0; listIndex < 2; ++listIndex )
    {
      NewRecord** hashTable = sp_newRecordHashTables[listIndex];
      
      if ( hashTable != 0 )
      {
	for ( unsigned int index = 0; index < bigPrime; ++index )
        {
	  NewRecord *p_thisRecord = hashTable[index];
        
	  while ( p_thisRecord != 0 ) {
	    if ( ! p_thisRecord->mb_deleted &&
		 p_thisRecord->msz_file != sz_unknownFile )
	      fprintf( stderr, "  %ld bytes\tallocated at %s:%d (%p)\n",
		       p_thisRecord->m_size,
		       p_thisRecord->msz_file,
		       p_thisRecord->m_line,
		       p_thisRecord->mp_data );
	    p_thisRecord = p_thisRecord->mp_next;
	  }
	}
      }
    }
    
    s_trackingOn = false;
    TracebackThisProcess();
  }
  else {
    fprintf( stderr, "\nAll bookends are intact.\n" );
  }
}
  

// A function called at exit to print out the undeleted records.

void
check_memory_tracker_records()
{
  cout << flush;
  cerr << flush;

  // If tracking has been turned off before we get here, the records
  // won't be consistent with reality.
  if ( ! s_trackingOn )
  {
    fprintf( stderr, "\nMemory tracking was turned off, so no list of undeleted records will be printed.\n" );
    return;
  }

  float avgUsage = 0.0;
  int   maxUsage = 0;
  int   numUsed  = 0;

  bool  foundUndeleted = false;

  for ( int listIndex = 0; listIndex < 2; ++listIndex )
  {
    NewRecord** hashTable = sp_newRecordHashTables[listIndex];

    if ( hashTable != 0 )
    {
      for ( unsigned int index = 0; index < bigPrime; ++index )
      {
        NewRecord *p_thisRecord = hashTable[index];
        NewRecord *p_deadRecord;
        
        int usage = 0;
        
        while ( p_thisRecord != 0 )
        {
          // We don't print out allocations of unknown origin, as we
          // can't ensure that all the static variables have been
          // destroyed at this point, and we don't want to stress out
          // the user needlessly. :)
          if ( ! p_thisRecord->mb_deleted &&
               p_thisRecord->msz_file != sz_unknownFile )
          {
            if ( ! foundUndeleted )
            {
              fprintf( stderr, "\nFound blocks that were allocated but never deleted:\n" );
              foundUndeleted = true;
            }
            fprintf( stderr, "  %ld bytes\tallocated at %s:%d (%p)\n",
                     p_thisRecord->m_size,
                     p_thisRecord->msz_file,
                     p_thisRecord->m_line,
		     p_thisRecord->mp_data );
	    if ( USE_BOOKENDS )
	      check_bookends( p_thisRecord );
          }
          p_deadRecord = p_thisRecord;
	  
          p_thisRecord = p_thisRecord->mp_next;
          ++usage;
        }
        
        if ( usage > 0 )
        {
          ++numUsed;
          avgUsage += ( (float)usage - avgUsage ) / (float)numUsed;
          if (maxUsage < usage )
            maxUsage = usage;
        }
      }
    }
  }

  if ( ! foundUndeleted ) {
    fprintf( stderr, "\nNo memory leaks detected.\n" );
  }

  //   if ( numUsed > 0 )
  //   {
  //     fprintf( stderr, "\n" );
  //     fprintf( stderr, "\nMemTracker stats:\n" );
  //     fprintf( stderr, "  number of bins used: %d\n", numUsed );
  //     fprintf( stderr, "  maximum size of used bin: %d\n", maxUsage );
  //     fprintf( stderr, "  average size of used bin: %f\n", avgUsage );
  //   }
}  

long
get_memory_usage()
{
  return memory_used;
}

#else

void
check_memory_tracker_records()
{
}

long
get_memory_usage()
{
  return 0;
}

#endif
