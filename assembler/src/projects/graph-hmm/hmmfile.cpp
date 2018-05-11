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

static int toESLAlphabet(Alphabet alph) {
    switch (alph) {
        case Alphabet::DNA:
            return eslDNA;
        case Alphabet::AMINO:
            return eslAMINO;
        default:
            return eslUNKNOWN;
    }

    return eslUNKNOWN;
}

static const char *toScoreSystem(ScoreSystem score, Alphabet alph) {
    switch (score) {
        case ScoreSystem::Default:
        default:
            return (alph == Alphabet::DNA ? "DNA1" : "BLOSUM62");
        case ScoreSystem::PAM30:
            return "PAM30";
        case ScoreSystem::PAM70:
            return "PAM70";
        case ScoreSystem::PAM120:
            return "PAM120";
        case ScoreSystem::PAM240:
            return "PAM240";
        case ScoreSystem::BLOSUM45:
            return "BLOSUM45";
        case ScoreSystem::BLOSUM50:
            return "BLOSUM50";
        case ScoreSystem::BLOSUM62:
            return "BLOSUM62";
        case ScoreSystem::BLOSUM80:
            return "BLOSUM80";
        case ScoreSystem::BLOSUM90:
            return "BLOSUM90";
    }
    return (alph == Alphabet::DNA ? "DNA1" : "BLOSUM62");
}

HMMSequenceBuilder::HMMSequenceBuilder(Alphabet alph, ScoreSystem score,
                                       double popen, double pext)
        : abc_(NULL, esl_alphabet_Destroy),
          bg_(NULL, p7_bg_Destroy),
          bld_(NULL, p7_builder_Destroy) {
    abc_.reset(esl_alphabet_Create(toESLAlphabet(alph)));
    bg_.reset(p7_bg_Create(abc_.get()));
    bld_.reset(p7_builder_Create(NULL, abc_.get()));

    p7_builder_LoadScoreSystem(bld_.get(), toScoreSystem(score, alph), popen, pext, bg_.get());
}

HMM HMMSequenceBuilder::from_string(const char *name, const char *seq, const char *desc) const {
    std::unique_ptr<ESL_SQ, void(*)(ESL_SQ*)> dbsq(esl_sq_CreateFrom(name, seq, desc, NULL, NULL),
                                                   esl_sq_Destroy);

    esl_sq_Digitize(abc_.get(), dbsq.get());

    return from_string(dbsq.get());
}

HMM HMMSequenceBuilder::from_string(ESL_SQ *dbsq) const {
    P7_HMM *hmm  = NULL;
    char errbuf[eslERRBUFSIZE];

    p7_SingleBuilder(bld_.get(), dbsq, bg_.get(), &hmm, NULL, NULL, NULL);

    if (p7_hmm_Validate(hmm, errbuf, 1e-5f) != eslOK) esl_fatal("HMM validation failed: %s\n", errbuf);

    ESL_ALPHABET *labc = esl_alphabet_Create(abc_->type);
    hmm->abc = labc;

    return { hmm, labc };
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
