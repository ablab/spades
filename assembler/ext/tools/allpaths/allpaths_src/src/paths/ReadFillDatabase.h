// Copyright (c) 2005 Broad Institute/Massachusetts Institute of Technology

#ifndef READ_FILL_DATABASE
#define READ_FILL_DATABASE

#include "system/System.h"
#include "Vec.h"
#include "paths/ReadFillRecord.h"

/// The file reads.{PREFIX}_fillrecords.k48 in the run_dir is a 
/// vec<ReadFillRecord> (one entry per read) created by CloseAllReadGaps,
/// CleanFills, or anything else that does filling.
///
/// Hand the filename to a ReadFillDatabase to read this information in 
/// an encapsulated way.  If you hand it a filename which does not exist, 
/// it will simulate the identity map.
///
/// A ReadFillDatabase also knows the left and right trim of each fill,
/// relative to the original read.  So the length of the trimmed read is
///   (orig read length in bases) - (left trim) - (right trim) - K + 1
/// Well, approximately, because of the flexible size of gaps.

class ReadFillDatabase {
public:

  ReadFillDatabase() 
    : mp_RFD( new identityRFD() ) {}

  // If the file exists, load data from there, otherwise get identity.
  ReadFillDatabase( const String& filename )
    : mp_RFD( NULL )    
  {
   FillsFromFile( filename );
  }
  
  ~ReadFillDatabase() { delete mp_RFD; }


  // This creates a vectorRFD (if filename exists) or an identityRFD (if not),
  // and loads trims, per-fill or per-read (also depending on what files exist).
  void FillsFromFile( const String& filename );

  // Pass these questions along to the pointed-to implementation

  int FirstFilling( int read ) const { return mp_RFD->FirstFilling(read); }
  int LastFilling( int read ) const { return mp_RFD->LastFilling(read); }
  int NumFillings( int read ) const { return mp_RFD->NumFillings(read); }
  int FillToRead( int read ) const { return mp_RFD->FillToRead(read); }
  bool Exploded( int read ) const { return mp_RFD->Exploded(read); }
  bool HitMergerLimit( int read ) const { return mp_RFD->HitMergerLimit(read); }
  bool IsDirty( int read ) const { return(Exploded(read) || HitMergerLimit(read)); }
  // Am I the trivial RFD?
  bool IsTrivialRFD() const { return mp_RFD->IsTrivialRFD(); }

  // Answer questions about trims here:
  int LeftTrimOfFill( int fill ) const {
    if( left_trim.empty() )
      return 0;
    else
      return left_trim[ trim_indexed_by_fill ? fill : FillToRead(fill) ];
  }
  int RightTrimOfFill( int fill ) const {
    if( right_trim.empty() )
      return 0;
    else
      return right_trim[ trim_indexed_by_fill ? fill : FillToRead(fill) ];
  }


private:
  class RFDimplementation;
  RFDimplementation* mp_RFD;
  bool trim_indexed_by_fill; // if false, trims are indexed by original read
  vec<int> left_trim;
  vec<int> right_trim;

private:
  // virtual base class, followed by two derived classes.
  class RFDimplementation {
  public:
    virtual int FirstFilling( int read ) const = 0;
    virtual int LastFilling( int read ) const = 0;
    virtual int NumFillings( int read ) const = 0;
    virtual int FillToRead( int read ) const = 0;
    virtual bool Exploded( int read ) const = 0;
    virtual bool HitMergerLimit( int read ) const = 0;
    virtual bool IsTrivialRFD() { return false; }

    RFDimplementation() {};
    virtual ~RFDimplementation() {};
  };

  // The identity ReadFillDatabase:
  class identityRFD : public RFDimplementation {
  public:
    identityRFD() : RFDimplementation() {};
    int FirstFilling( int read ) const { return read; }
    int LastFilling( int read ) const { return read; }
    int NumFillings( int read ) const { return 1; }
    int FillToRead( int read ) const { return read; }
    bool Exploded( int read ) const { return false; }
    bool HitMergerLimit( int read ) const { return false; }
    bool IsTrivialRFD() { return true; }
  };

  // An interface to a vec<ReadFillRecord>
  class vectorRFD : public RFDimplementation {
  public:
    vectorRFD() 
      : RFDimplementation(), 
	mp_fills( NULL ),
	I_own_my_readfills( false ),
	cached_read( 0 )
    {}

    vectorRFD( const String& filename )
      : RFDimplementation(), 
	I_own_my_readfills( false ), // changed to true by FillsFromFile
	cached_read( 0 )
    { FillsFromFile(filename); }

    ~vectorRFD() { if( I_own_my_readfills ) delete mp_fills; }

    void FillsFromFile( String filename ) {
      if( I_own_my_readfills ) delete mp_fills;
      vec<ReadFillRecord>* p_fills = new vec<ReadFillRecord>;
      BinaryRead2( filename, *p_fills );
      mp_fills = p_fills;
      I_own_my_readfills = true;
      cached_read = 0;
    }

    // Queries we pass along to the appropriate ReadFillRecord.
    // These are all constant-time.

    int FirstFilling( int read ) const 
    { return (*mp_fills)[read].FirstFilling(); }

    int LastFilling( int read ) const 
    { return (*mp_fills)[read].LastFilling(); }

    int NumFillings( int read ) const 
    { return (*mp_fills)[read].NumFillings(); }

    bool Exploded( int read ) const 
    { return (*mp_fills)[read].Exploded(); }

    bool HitMergerLimit( int read ) const 
    { return (*mp_fills)[read].HitMergerLimit(); }

    // Which read did the ith fill come from?
    // Requires a binary search of the vector, so logarithmic time.
    // But we cache the most recent result.  This helps when you want to get 
    // the trims for all the fillings of a read whose length never changed.

    int FillToRead( int fill ) const;

  private:
    const vec<ReadFillRecord>* mp_fills;
    bool I_own_my_readfills;  // hold-over: now this is always true
    mutable int cached_read;  // cache for FillToRead binray search
  };  // end of class vectorRFD : public RFDimplementation

};  // end of class ReadFillDatabase

#endif
