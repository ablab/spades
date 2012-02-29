///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#ifndef SHARED_MEM_H
#define SHARED_MEM_H

#include <sys/mman.h>

#include "CoreTools.h"

class SharedMem {

     public:

     SharedMem( ) : owner_(False) { }

     void Initialize( const String& sharename, const int len = 0, 
          Bool owner = False ) 
     {    sharename_ = sharename;
          len_ = len;
          owner_ = owner;
          int fd = shm_open( sharename.c_str( ), 
                ( owner ? O_RDWR ^ O_CREAT : O_RDWR ), 0664 );
          ForceAssertGe( fd, 0 );
          if ( len > 0 && owner )
          {    int fstatus = ftruncate( fd, len );
               ForceAssertEq( fstatus, 0 );    }
          void* pa = mmap( 0, len, PROT_READ ^ PROT_WRITE, MAP_SHARED, fd, 0 );
          Close(fd);
          loc_ = reinterpret_cast<unsigned char*>(pa);    }

     SharedMem( const String& sharename, const int len = 0, 
          const Bool owner = False ) 
     {    Initialize( sharename, len, owner );    }

     void resize( const int len )
     {    int fd = shm_open( sharename_.c_str( ), O_RDWR, 0664 );
          ForceAssertGe( fd, 0 );
          int fstatus = ftruncate( fd, len );
          ForceAssertEq( fstatus, 0 );
          void* pa = mmap( 0, len, PROT_READ ^ PROT_WRITE, MAP_SHARED, fd, 0 );
          loc_ = reinterpret_cast<unsigned char*>(pa);
          Close(fd);
          len_ = len;    }

     void SetLen( const int len )
     {    len_ = len;    }

     unsigned char* Loc( )
     {    return loc_;    }

     void SetLoc( unsigned char* loc )
     {    loc_ = loc;    }

     ~SharedMem( )
     {    if (owner_) shm_unlink( sharename_.c_str( ) );    }

     unsigned char& operator[ ]( const int i )
     {    return loc_[i];    }

     String ShareName( ) const { return sharename_; }
     int Len( ) const { return len_; }
     volatile const unsigned char* Loc( ) const { return loc_; }

     private:

     String sharename_;
     int len_;
     unsigned char* loc_;
     Bool owner_;

};

#endif
