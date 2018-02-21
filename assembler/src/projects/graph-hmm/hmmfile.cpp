#include "hmmfile.hpp"

extern "C" {
#include "hmmer.h"
};

namespace hmmer {
HMMFile::HMMFile(const std::string &hmmfile)
        : hmmfile_(NULL, p7_hmmfile_Close) {
    int status   = eslOK;
    P7_HMMFILE *hfp  = NULL;

    status = p7_hmmfile_OpenE(hmmfile.c_str(), nullptr, &hfp, NULL);
    if (status == eslOK)
        hmmfile_.reset(hfp);

    //if      (status == eslENOTFOUND) p7_Fail("File existence/permissions problem in trying to open HMM file %s.\n%s\n", cfg.hmmfile.c_str(), errbuf);
    //else if (status == eslEFORMAT)   p7_Fail("File format problem in trying to open HMM file %s.\n%s\n",                cfg.hmmfile.c_str(), errbuf);
    //maelse if (status != eslOK)        p7_Fail("Unexpected error %d in opening HMM file %s.\n%s\n",               status, cfg.hmmfile.c_str(), errbuf);
}

HMM::HMM(P7_HMM *hmm, ESL_ALPHABET *abc) 
        : hmm_(hmm, p7_hmm_Destroy),
          abc_(abc, esl_alphabet_Destroy)
{}

static ReadErrc to_ReadErrc(int err) {
    switch (err) {
        case eslOK:
            return ReadErrc::NoError;
        case eslEOD:
            return ReadErrc::EOD;
        case eslEFORMAT:
            return ReadErrc::BadFormat;
        case eslEINCOMPAT:
            return ReadErrc::ABCIncompat;
        case eslEOF:
            return ReadErrc::EndOfFile;
    };

    return ReadErrc::Unknown;
}

ErrorOr<HMM> HMMFile::read() {
    ESL_ALPHABET *abc = NULL;
    P7_HMM *rhmm = NULL;
    
    int res = p7_hmmfile_Read(hmmfile_.get(), &abc, &rhmm);
    if (res != eslOK)
        return to_ReadErrc(res);

    return HMM(rhmm, abc);
}

}  // namespace hmmer

namespace {
struct HMMReadErrCategory : std::error_category {
  const char* name() const noexcept override;
  std::string message(int ev) const override;
};

const char* HMMReadErrCategory::name() const noexcept {
  return "HMM read";
}

std::string HMMReadErrCategory::message(int ev) const {
    switch (static_cast<hmmer::ReadErrc>(ev)) {
        case hmmer::ReadErrc::EOD:
            return "read failed, HMM file may be truncated?";
        case hmmer::ReadErrc::BadFormat:
            return "bad file format in HMM file";
        case hmmer::ReadErrc::ABCIncompat:
            return "incompatible HMM alphabet";
        case hmmer::ReadErrc::EndOfFile:
            
        default:
            return "(unknown error)";
    }
}

const HMMReadErrCategory theHMMReadErrCategory {};

}  // namespace

namespace hmmer {
std::error_code make_error_code(ReadErrc e) {
  return {static_cast<int>(e), theHMMReadErrCategory};
}
}




