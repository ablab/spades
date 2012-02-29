///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#include "feudal/FeudalTools.h"
#include "feudal/MasterVec.h"
#include "feudal/FeudalControlBlock.h"
#include "String.h"
#include "Vec.h"

namespace
{
class MVCopier
{
public:
    MVCopier( String const& filename )
    : mIS(filename.c_str()), mFCB(0,0,0,0,0)
    { struct stat sb;
      if ( stat(filename.c_str(),&sb) )
      { FatalErr("Can't stat file " << filename); }
      if ( !mIS.read(reinterpret_cast<char*>(&mFCB),sizeof(mFCB)) )
      { FatalErr("Can't read FCB from " << filename); }
      if ( !mFCB.isValid(filename.c_str(),sb.st_size,true) )
      { FatalErr("FCB invalid for " << filename); }
      if ( mFCB.isCompressed() )
      { FatalErr("Can't merge compressed mastervec files."); }
      mFileSize = sb.st_size; }

    size_t getNElements() const { return mFCB.getNElements(); }
    size_t getVarDataLen() const { return mFCB.getVarDataLen(); }
    size_t getFixedDataLen() const { return mFileSize - mFCB.getFixedOffset(); }
    FeudalControlBlock const& getFCB() const { return mFCB; }

    void copy( ostream& os, size_t len )
    { char buf[8192];
      while ( len )
      { size_t nnn = std::min(len,sizeof(buf));
        if ( !mIS.read(buf,nnn) )
        { FatalErr("Can't read var-data input."); }
        if ( !os.write(buf,nnn) )
        { FatalErr("Can't write var-data."); }
        len -= nnn; } }

    size_t readOffset()
    { size_t result;
      if ( !mIS.read(reinterpret_cast<char*>(&result),sizeof(result)) )
      { FatalErr("Can't read offset."); }
      return result; }

private:
    MVCopier( MVCopier const& );
    MVCopier& operator=( MVCopier const& );

    ifstream mIS;
    FeudalControlBlock mFCB;
    size_t mFileSize;
};

} // end of anonymous namespace


bool IsGoodFeudalFile( const String& filename, bool verbose )
{
    return FeudalControlBlock::isGoodFeudalFile(filename.c_str(),verbose);
}

void RemoveMastervecFiles(const String& filename) {
  Remove( filename );
  Remove( filename + ".gz" );
  Remove( filename + "..offsets" );
  Remove( filename + "..static" );
}

void MergeMastervecs( const String& dstFile, std::vector<String> const& inputs )
{
    if ( inputs.size() < 2 )
    { FatalErr("Can't merge fewer than 2 feudal files."); }

    std::vector<MVCopier*> copiers;
    copiers.reserve(inputs.size());

    size_t totElements = 0;
    size_t totVarData = 0;
    std::vector<String>::const_iterator namesEnd(inputs.end());
    for ( std::vector<String>::const_iterator namesItr(inputs.begin());
            namesItr != namesEnd; ++namesItr )
    {
        MVCopier* pCopier = new MVCopier(*namesItr);
        totElements += pCopier->getNElements();
        totVarData += pCopier->getVarDataLen();
        copiers.push_back(pCopier);
    }

    ofstream os(dstFile.c_str());

    FeudalControlBlock const& hdr1 = copiers.front()->getFCB();
    FeudalControlBlock dstHdr(totElements,totVarData,
                              hdr1.getSizeofFixed(),
                              hdr1.getSizeofX(),
                              hdr1.getSizeofA(),
                              1, false);

    if ( !os.write(reinterpret_cast<char*>(&dstHdr),sizeof(dstHdr)) )
    { FatalErr("Can't write FCB to " << dstFile); }

    // move the var data
    std::vector<MVCopier*>::iterator end(copiers.end());
    for ( std::vector<MVCopier*>::iterator itr(copiers.begin()); itr != end; ++itr )
        (*itr)->copy(os,(*itr)->getVarDataLen());

    // move the var-data offset tables

    size_t lastOffset = sizeof(dstHdr);
    if ( !os.write(reinterpret_cast<char*>(&lastOffset),sizeof(lastOffset)) )
    { FatalErr("Can't write first offset-table element."); }

    for ( std::vector<MVCopier*>::iterator itr(copiers.begin()); itr != end; ++itr )
    {
        long offsetDiff = lastOffset - (*itr)->readOffset();
        size_t nnn = (*itr)->getNElements();
        while ( nnn-- )
        {
            lastOffset = (*itr)->readOffset() + offsetDiff;
            if ( !os.write(reinterpret_cast<char*>(&lastOffset),sizeof(lastOffset)) )
            { FatalErr("Can't write offset-table element."); }
        }
    }

    // move the fixed-length data
    for ( std::vector<MVCopier*>::iterator itr(copiers.begin()); itr != end; ++itr )
        (*itr)->copy(os,(*itr)->getFixedDataLen());

    os.close();

    for ( std::vector<MVCopier*>::iterator itr(copiers.begin()); itr != end; ++itr )
        delete *itr;
}

void MergeMastervecs(const String& file1, const String& file2, const String& dest_file)
{
    std::vector<String> inputs;
    inputs.push_back(file1);
    inputs.push_back(file2);
    MergeMastervecs(dest_file, inputs);
}


size_t MastervecFileObjectCount(const String& filename)
{
    FeudalControlBlock fcb(filename.c_str(), false);
    return fcb.getNElements();
}

size_t MastervecFileRawCount( const String& filename, size_t dataSize )
{
  if ( !IsRegularFile(filename) && !IsRegularFile( filename + ".gz" ) )
    FatalErr( "Neither " << filename << " nor " << filename << ".gz exists." );
  if ( IsRegularFile(filename) && IsRegularFile( filename + ".gz" ) )
    FatalErr( "Both " << filename << " and " << filename << ".gz exist."
	      << "  This confuses me." );

  if ( IsRegularFile( filename + ".gz" ) )
    System( "gzip -d " + filename + ".gz" );

  FeudalControlBlock fcb(filename.c_str(),false);

  //get the data size from the file header, if we can.
  if ( 0 == dataSize ) dataSize = fcb.getSizeofA();
  Assert(dataSize);

  if ( fcb.getNFiles() != 1 && fcb.getNFiles() != 3 )
    FatalErr( "ReadAll: It would appear that " << filename << "\n"
	      << "is not in proper mastervec format.\n"
	      << "The value of nfiles is " << fcb.getNFiles()
	      << ", whereas it should be 1 or 3,\nand for ReadAll, it should "
	      << "be 1.\nThere is probably something very wrong with the file." );
  if ( fcb.getNFiles() != 1 )
    FatalErr( "ReadAll: " << filename << " was expected to be in "
	      << "mastervec single-file format, but it has "
	      << "control.nfiles() = 3." );

  return fcb.getVarDataLen()/dataSize;
}

