#include "ErrorOr.hpp"

#include <memory>
#include <string>
#include <system_error>

#include "hmmer_fwd.h"

namespace hmmer {

class HMMFile;

class HMM {
  public:
    HMM(P7_HMM *hmm = NULL,
        ESL_ALPHABET *abc = NULL);

    P7_HMM *get() const { return hmm_.get(); }
    ESL_ALPHABET *abc() const { return abc_.get(); }
    
  private:
    friend class HMMFile;
    std::unique_ptr<P7_HMM, void(*)(P7_HMM*)> hmm_;
    std::unique_ptr<ESL_ALPHABET, void(*)(ESL_ALPHABET*)> abc_;
};

enum class OpenErrc {
    NoError = 0,
    NotFound,
    BadFormat,
    UnknownError
};

enum class ReadErrc {
    NoError = 0,
    EOD,
    BadFormat,
    ABCIncompat,
    EndOfFile,
    Unknown
};

enum class OpenErrc;
std::error_code make_error_code(OpenErrc);
enum class ReadErrc;
std::error_code make_error_code(ReadErrc);

class HMMFile {
  public:
    HMMFile(const std::string &name);

    ErrorOr<HMM> read();
    bool valid() const { return (bool)hmmfile_; }
    P7_HMMFILE *get() const { return hmmfile_.get(); }
    
  private:
    std::unique_ptr<P7_HMMFILE, void(*)(P7_HMMFILE*)> hmmfile_;
};

enum class Alphabet {
    DNA,
    AMINO
};

enum class ScoreSystem {
    Default,
    PAM30,
    PAM70,
    PAM120,
    PAM240,
    BLOSUM45,
    BLOSUM50,
    BLOSUM62,
    BLOSUM80,
    BLOSUM90
};

class HMMSequenceBuilder {
  public:
    HMMSequenceBuilder(Alphabet alph, ScoreSystem score,
                       double popen = 0.1, double pext = 0.4);

    HMM from_string(const char *name, const char *seq, const char *desc) const;
    HMM from_string(ESL_SQ *dbsq) const;
    
  private:
    std::unique_ptr<ESL_ALPHABET, void(*)(ESL_ALPHABET*)> abc_;
    std::unique_ptr<P7_BG, void(*)(P7_BG*)> bg_;
    std::unique_ptr<P7_BUILDER, void(*)(P7_BUILDER*)> bld_;
};


}  // namespace hmmer

namespace std {
template <>
struct is_error_code_enum<hmmer::OpenErrc> : true_type {};
template <>
struct is_error_code_enum<hmmer::ReadErrc> : true_type {};
}

