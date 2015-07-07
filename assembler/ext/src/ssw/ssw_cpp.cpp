#include "ssw/ssw_cpp.h"

#include <sstream>

extern "C" {
#include "ssw/ssw.h"
}

namespace {

static int8_t kBaseTranslation[128] = {
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
  //   A     C            G
    4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4, 
  //             T
    4, 4, 4, 4,  3, 0, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
  //   a     c            g
    4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4, 
  //             t
    4, 4, 4, 4,  3, 0, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4 
};

void BuildSwScoreMatrix(const uint8_t& match_score, 
                        const uint8_t& mismatch_penalty,
			int8_t* matrix) {

  // The score matrix looks like
  //                 // A,  C,  G,  T,  N
  //  score_matrix_ = { 2, -2, -2, -2,  0, // A
  //                   -2,  2, -2, -2,  0, // C
  //                   -2, -2,  2, -2,  0, // G
  //                   -2, -2, -2,  2,  0, // T
  //                    0,  0,  0,  0,  0};// N

  int id = 0;
  for (int i = 0; i < 4; ++i) {
    for (int j = 0; j < 4; ++j) {
      matrix[id] = ((i == j) ? match_score : static_cast<int8_t>(-mismatch_penalty));
      ++id;
    }
    matrix[id] = 0;
    ++id;
  }

  for (int i = 0; i < 5; ++i)
    matrix[id++] = 0;
    
}

void ConvertAlignment(const s_align& s_al, 
                      const int& query_len, 
                      StripedSmithWaterman::Alignment* al) {
  al->sw_score           = s_al.score1;
  al->sw_score_next_best = s_al.score2;
  al->ref_begin          = s_al.ref_begin1;
  al->ref_end            = s_al.ref_end1;
  al->query_begin        = s_al.read_begin1;
  al->query_end          = s_al.read_end1;
  al->ref_end_next_best  = s_al.ref_end2;

  al->cigar.clear();
  al->cigar_string.clear();
  
  if (s_al.cigarLen > 0) {
    std::ostringstream cigar_string;
    if (al->query_begin > 0) {
      uint32_t cigar = (al->query_begin << 4) | 0x0004;
      al->cigar.push_back(cigar);
      cigar_string << al->query_begin << 'S';
    }

    for (int i = 0; i < s_al.cigarLen; ++i) {
      al->cigar.push_back(s_al.cigar[i]);
      cigar_string << (s_al.cigar[i] >> 4);
      uint8_t op = s_al.cigar[i] & 0x000f;
      switch(op) {
        case 0: cigar_string << 'M'; break;
        case 1: cigar_string << 'I'; break;
        case 2: cigar_string << 'D'; break;
      }
    }

    int end = query_len - al->query_end - 1;
    if (end > 0) {
      uint32_t cigar = (end << 4) | 0x0004;
      al->cigar.push_back(cigar);
      cigar_string << end << 'S';
    }

    al->cigar_string = cigar_string.str();
  } // end if
}

int CalculateNumberMismatch(
    const StripedSmithWaterman::Alignment& al,
    const int8_t* matrix,
    int8_t const *ref,
    int8_t const *query) {
  
  ref   += al.ref_begin;
  query += al.query_begin;
  int mismatch_length = 0;
  for (unsigned int i = 0; i < al.cigar.size(); ++i) {
    int32_t op = al.cigar[i] & 0x0000000f;
    int32_t length = (al.cigar[i] >> 4) & 0x0fffffff;
    if (op == 0) { // M
      for (int j = 0; j < length; ++j) {
        if (matrix[*ref] != matrix[*query]) ++mismatch_length;
	++ref;
	++query;
      }
    } else if (op == 1) { // I
      query += length;
      mismatch_length += length;
    } else if (op == 2) { // D
      ref += length;
      mismatch_length += length;
    }
  }

  return mismatch_length;
}

void SetFlag(const StripedSmithWaterman::Filter& filter, uint8_t* flag) {
  if (filter.report_begin_position) *flag |= 0x08;
  if (filter.report_cigar) *flag |= 0x0f;
}

} // namespace



namespace StripedSmithWaterman {

Aligner::Aligner(void)
    : score_matrix_(NULL)
    , score_matrix_size_(5)
    , translation_matrix_(NULL)
    , default_matrix_(false)
    , matrix_built_(false)
    , match_score_(2)
    , mismatch_penalty_(2)
    , gap_opening_penalty_(3)
    , gap_extending_penalty_(1)
    , translated_reference_(NULL)
    , reference_length_(0)
{
  BuildDefaultMatrix();
}

Aligner::Aligner(
    const uint8_t& match_score,
    const uint8_t& mismatch_penalty,
    const uint8_t& gap_opening_penalty,
    const uint8_t& gap_extending_penalty)

    : score_matrix_(NULL)
    , score_matrix_size_(5)
    , translation_matrix_(NULL)
    , default_matrix_(false)
    , matrix_built_(false)
    , match_score_(match_score)
    , mismatch_penalty_(mismatch_penalty)
    , gap_opening_penalty_(gap_opening_penalty)
    , gap_extending_penalty_(gap_extending_penalty)
    , translated_reference_(NULL)
    , reference_length_(0)
{
  BuildDefaultMatrix();
}

Aligner::Aligner(const int8_t* score_matrix,
                 const int&    score_matrix_size,
	         const int8_t* translation_matrix,
		 const int&    translation_matrix_size)
    
    : score_matrix_(NULL)
    , score_matrix_size_(score_matrix_size)
    , translation_matrix_(NULL)
    , default_matrix_(true)
    , matrix_built_(false)
    , match_score_(2)
    , mismatch_penalty_(2)
    , gap_opening_penalty_(3)
    , gap_extending_penalty_(1)
    , translated_reference_(NULL)
    , reference_length_(0)
{
  score_matrix_ = new int8_t[score_matrix_size_ * score_matrix_size_];
  memcpy(score_matrix_, score_matrix, sizeof(int8_t) * score_matrix_size_ * score_matrix_size_);
  translation_matrix_ = new int8_t[translation_matrix_size];
  memcpy(translation_matrix_, translation_matrix, sizeof(int8_t) * translation_matrix_size);
  matrix_built_ = true;
}


Aligner::~Aligner(void){
  Clear();
}

int Aligner::SetReferenceSequence(const char* seq, const int& length) {
  
  int len = 0;
  if (matrix_built_) {
    // calculate the valid length
    //int calculated_ref_length = static_cast<int>(strlen(seq));
    //int valid_length = (calculated_ref_length > length) 
    //                   ? length : calculated_ref_length;
    int valid_length = length;
    // delete the current buffer
    CleanReferenceSequence();
    // allocate a new buffer
    translated_reference_ = new int8_t[valid_length];
  
    len = TranslateBase(seq, valid_length, translated_reference_);
  } else {
    // nothing
  }

  reference_length_ = len;
  return len;


}

int Aligner::TranslateBase(const char* bases, const int& length, 
    int8_t* translated) const {

  char* ptr = (char*)bases;
  int len = 0;
  for (int i = 0; i < length; ++i) {
    translated[i] = translation_matrix_[(int) *ptr];
    ++ptr;
    ++len;
  }

  return len;
}


bool Aligner::Align(const char* query, const Filter& filter, 
                    Alignment* alignment) const
{
  if (!matrix_built_) return false;
  if (reference_length_ == 0) return false;

  int query_len = strlen(query);
  if (query_len == 0) return false;
  int8_t* translated_query = new int8_t[query_len];
  TranslateBase(query, query_len, translated_query);

  const int8_t score_size = 2;
  s_profile* profile = ssw_init(translated_query, query_len, score_matrix_, 
                                score_matrix_size_, score_size);

  uint8_t flag = 0;
  SetFlag(filter, &flag);
  s_align* s_al = ssw_align(profile, translated_reference_, reference_length_,
                                 static_cast<int>(gap_opening_penalty_), 
				 static_cast<int>(gap_extending_penalty_),
				 flag, filter.score_filter, filter.distance_filter, query_len);
  
  alignment->Clear();
  ConvertAlignment(*s_al, query_len, alignment);
  alignment->mismatches = CalculateNumberMismatch(*alignment, score_matrix_, translated_reference_, translated_query);


  // Free memory
  if (query_len > 1) delete [] translated_query;
  else delete translated_query;
  align_destroy(s_al);
  init_destroy(profile);

  return true;
}


bool Aligner::Align(const char* query, const char* ref, const int& ref_len,
                    const Filter& filter, Alignment* alignment) const
{
  if (!matrix_built_) return false;
  
  int query_len = strlen(query);
  if (query_len == 0) return false;
  int8_t* translated_query = new int8_t[query_len];
  TranslateBase(query, query_len, translated_query);

  // calculate the valid length
  //int calculated_ref_length = static_cast<int>(strlen(ref));
  //int valid_ref_len = (calculated_ref_length > ref_len) 
  //                    ? ref_len : calculated_ref_length;
  int valid_ref_len = ref_len;
  int8_t* translated_ref = new int8_t[valid_ref_len];
  TranslateBase(ref, valid_ref_len, translated_ref);


  const int8_t score_size = 2;
  s_profile* profile = ssw_init(translated_query, query_len, score_matrix_, 
                                score_matrix_size_, score_size);

  uint8_t flag = 0;
  SetFlag(filter, &flag);
  s_align* s_al = ssw_align(profile, translated_ref, valid_ref_len,
                                 static_cast<int>(gap_opening_penalty_), 
				 static_cast<int>(gap_extending_penalty_),
				 flag, filter.score_filter, filter.distance_filter, query_len);
  
  alignment->Clear();
  ConvertAlignment(*s_al, query_len, alignment);
  alignment->mismatches = CalculateNumberMismatch(*alignment, score_matrix_, translated_ref, translated_query);

  // Free memory
  if (query_len > 1) delete [] translated_query;
  else delete translated_query;
  if (valid_ref_len > 1) delete [] translated_ref;
  else delete translated_ref;
  align_destroy(s_al);
  init_destroy(profile);

  return true;
}

void Aligner::Clear(void) {
  if (score_matrix_) delete [] score_matrix_;
  score_matrix_ = NULL;

  if (!default_matrix_ && translation_matrix_) 
    delete [] translation_matrix_;
  translation_matrix_ = NULL;

  CleanReferenceSequence();

  default_matrix_ = false;
  matrix_built_   = false;
}

void Aligner::SetAllDefault(void) {
  score_matrix_size_     = 5;
  default_matrix_        = false;
  matrix_built_          = false;
  match_score_           = 2;
  mismatch_penalty_      = 2;
  gap_opening_penalty_   = 3;
  gap_extending_penalty_ = 1;
  reference_length_      = 0;
}

bool Aligner::ReBuild(void) {
  if (matrix_built_) return false;

  SetAllDefault();
  BuildDefaultMatrix();

  return true;
}

bool Aligner::ReBuild(
    const uint8_t& match_score,
    const uint8_t& mismatch_penalty,
    const uint8_t& gap_opening_penalty,
    const uint8_t& gap_extending_penalty) {
  if (matrix_built_) return false;

  SetAllDefault();

  match_score_           = match_score;
  mismatch_penalty_      = mismatch_penalty;
  gap_opening_penalty_   = gap_opening_penalty;
  gap_extending_penalty_ = gap_extending_penalty;

  BuildDefaultMatrix();

  return true;
}

bool Aligner::ReBuild(
    const int8_t* score_matrix,
    const int&    score_matrix_size,
    const int8_t* translation_matrix,
    const int&    translation_matrix_size) {

  score_matrix_ = new int8_t[score_matrix_size_ * score_matrix_size_];
  memcpy(score_matrix_, score_matrix, sizeof(int8_t) * score_matrix_size_ * score_matrix_size_);
  translation_matrix_ = new int8_t[translation_matrix_size];
  memcpy(translation_matrix_, translation_matrix, sizeof(int8_t) * translation_matrix_size);
  matrix_built_ = true;

  return true;
}

void Aligner::BuildDefaultMatrix(void) {
  score_matrix_ = new int8_t[score_matrix_size_ * score_matrix_size_];
  BuildSwScoreMatrix(match_score_, mismatch_penalty_, score_matrix_);
  translation_matrix_ = kBaseTranslation;
  matrix_built_   = true;
  default_matrix_ = true;
}
} // namespace StripedSmithWaterman
