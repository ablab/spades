///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2009) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////


#include "kmers/KmerParcels.h"




double NumSigmaBinomial(size_t n, size_t m, double p)
{
  const double mu = p * n;
  const double sig = sqrt((1 - p) * mu);
  return fabs(m - mu) / sig;
}  


inline String ANSI256RGBColorID(int red, int green, int blue)
{  return "8;5;" + ToString(16 + blue + 6 * (green + 6 * red)); }

inline String ANSISetBackground(String color) { return "\x1b[4" + color + "m"; }
inline String ANSISetForeground(String color) { return "\x1b[3" + color + "m"; }
inline String ANSISetDefault() { return "\x1b[0m"; }

static const char hiero[] = { '^', '(', '-', '.' };


void PrintKmerBaseVec(const BaseVec & bv, 
                      const size_t color, ostream& out)
{
  String colors[4];
  colors[0] = ANSI256RGBColorID(5, 0, 0);
  colors[1] = ANSI256RGBColorID(0, 0, 5);
  colors[2] = ANSI256RGBColorID(5, 5, 0);
  colors[3] = ANSI256RGBColorID(0, 5, 0);


  for (size_t ib = 0; ib != bv.size(); ib++) {
    const uint8_t b = bv[ib];
    const char base = (color == 3) ? hiero[b] : as_base(b);

    out << ANSISetForeground(colors[b]) + base;
  }
     
  out << ANSISetDefault();
}



void PrintKmerBaseVec(const BaseVec & bv, const QualNibbleVec & qv, 
                      const size_t color, ostream& out)
{
  for (size_t ib = 0; ib != bv.size(); ib++) {
    const uint8_t b = bv[ib];
      
    const char base = (color == 3) ? hiero[b] : as_base(b);

    size_t qual_color = Min(qv[ib] / 8, 5);
    String color_str = ANSI256RGBColorID(5 - qual_color, qual_color, 0);

    out << ANSISetForeground(color_str) + base;
  }
     
  out << ANSISetDefault();
}



// ---------------------------------------------
// KmerParcels    Store
// ---------------------------------------------


vec<size_t> KmerParcelsStore::GetParcelsIDsSizeSorted() const
{
  const size_t n_parcels = GetNumParcels();
  ForceAssertGt(n_parcels, 0u);
  
  vec<size_t> parcel_IDs(n_parcels);
  vec<size_t> n_batches(n_parcels);

  for (size_t parcel_ID = 0; parcel_ID != n_parcels; parcel_ID++) {
    KmerParcelReader parcel(*this, parcel_ID);
    n_batches[parcel_ID] = parcel.GetNumKmerBatches();
    parcel_IDs[parcel_ID] = parcel_ID;
  }
  ReverseSortSync(n_batches, parcel_IDs);
  
  return parcel_IDs;
}






size_t KmerParcelsStore::GetTotalNumKmerBatches() const
{
  const size_t n_parcels = GetNumParcels();
  ForceAssertGt(n_parcels, 0u);

  size_t n_batches = 0;

  for (size_t parcel_ID = 0; parcel_ID != n_parcels; parcel_ID++) {
    KmerParcelReader parcel(*this, parcel_ID);
    n_batches += parcel.GetNumKmerBatches();
  }
  
  return n_batches;
}





































