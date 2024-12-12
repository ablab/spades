//***************************************************************************
//* Copyright (c) 2024 SPAdes team
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "single_read.hpp"

#include "utils/verify.hpp"
#include "io/reads/parser.hpp"
#include "sequence/quality.hpp"
#include "sequence/nucl.hpp"

#include "kseq/kseq.h"

#include <zstd.h>
#include <string>
#include <cstdio>

namespace io {

namespace fastafastqzstd {

struct zstdFile_s {
    FILE *fin;
    ZSTD_DCtx *dctx;
    void *bufIn;
    size_t bufInSize;
    void *bufOut;
    size_t bufOutSize;

    ZSTD_inBuffer input;
    ZSTD_outBuffer output;
    size_t readPos;
};

typedef struct zstdFile_s *zstdFile;

size_t zstdread(zstdFile_s *state, void *buf, size_t len) {
    // See if there is something in the output buffer that we can reuse
    if (state->readPos >= state->output.pos) { // nope
        size_t no_progress = 0;
        do {
            // See if there is something in the input buffer
            if (state->input.pos >= state->input.size) { // nope
                // Read another input frame
                size_t read = fread(state->bufIn, 1, state->bufInSize, state->fin);
                if (read == 0)
                    return 0; // Nothing left anywhere
                state->input = { state->bufIn, read, 0 };
            }

            // There is something in input buffer, go ahead and decompress it
            state->output = { state->bufOut, state->bufOutSize, 0 };
            state->readPos = 0;
            size_t ret = ZSTD_decompressStream(state->dctx, &state->output, &state->input);
            VERIFY_MSG(!ZSTD_isError(ret),
                       "zstd decompression error, code: " << ret <<
                       ", decription: " << ZSTD_getErrorName(ret));
            VERIFY_MSG(no_progress++ < 16, "zlib decompression problem, no progress after " << no_progress << " decompression calls");
        } while (state->output.pos == 0);
    }

    len = std::min(len, state->output.pos - state->readPos);
    VERIFY(len > 0);
    memcpy(buf, (uint8_t*)state->output.dst + state->readPos, len);
    state->readPos += len;
    return len;
}

// STEP 1: declare the type of file handler and the read() function
// Silence bogus gcc warnings
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wconversion"
// STEP 1: declare the type of file handler and the read() function
KSEQ_INIT(zstdFile, zstdread)
#pragma GCC diagnostic pop
}

class FastaFastqZstdParser: public Parser {
public:
    /*
     * Default constructor.
     *
     * @param filename The name of the file to be opened.
     * @param offset The offset of the read quality.
     */
    FastaFastqZstdParser(const std::filesystem::path& filename,
                         FileReadFlags flags = FileReadFlags())
            : Parser(filename, flags), fp_(), seq_(NULL) {
        open();
    }

    /*
     * Default destructor.
     */
    /* virtual */
    ~FastaFastqZstdParser() {
        close();
    }

    /*
     * Read SingleRead from stream.
     *
     * @param read The SingleRead that will store read data.
     *
     * @return Reference to this stream.
     */
    /* virtual */
    FastaFastqZstdParser& operator>>(SingleRead& read) {
        if (!is_open_ || eof_)
            return *this;

        if (seq_->qual.s && flags_.use_name && flags_.use_quality) {
            read = SingleRead(seq_->name.s ? seq_->name.s : "",
                              seq_->comment.s ? seq_->comment.s : "",
                              seq_->seq.s, seq_->qual.s,
                              flags_.offset,
                              0, 0, flags_.validate);
        } else if (flags_.use_name && seq_->name.s) {
            read = SingleRead(seq_->name.s,
                              seq_->comment.s ? seq_->comment.s : "",
                              seq_->seq.s,
                              0, 0, flags_.validate);
        } else if (flags_.use_comment && seq_->comment.s) {
            read = SingleRead("", seq_->comment.s,
                              seq_->seq.s,
                              0, 0, flags_.validate);
        } else
            read = SingleRead(seq_->seq.s,
                              0, 0, flags_.validate);

        ReadAhead();
        return *this;
    }

    /*
     * Close the stream.
     */
    /* virtual */
    void close() {
        if (!is_open_)
            return;

        // STEP 5: destroy seq
        fastafastqzstd::kseq_destroy(seq_);
        // STEP 6: close the file handler
        fclose(fp_.fin);
        ZSTD_freeDStream(fp_.dctx);
        free(fp_.bufIn);
        free(fp_.bufOut);
        is_open_ = false;
        eof_ = true;
    }

private:
    /*
     * @variable File that is associated with gzipped data file.
     */
    fastafastqzstd::zstdFile_s fp_;
    /*
     * @variable Data element that stores last SingleRead got from
     * stream.
     */
    fastafastqzstd::kseq_t* seq_;

    /*
     * Open a stream.
     */
    /* virtual */
    void open() {
        // STEP 2: open the file handler
        fp_.fin = fopen(filename_.c_str(), "rb");
        if (!fp_.fin) {
            is_open_ = false;
            return;
        }
        
        // Initialize zstd context and all buffers
        fp_.bufInSize = ZSTD_DStreamInSize();
        fp_.bufOutSize = ZSTD_DStreamOutSize();
        fp_.bufIn = malloc(fp_.bufInSize);
        fp_.bufOut = malloc(fp_.bufOutSize);
        fp_.dctx = ZSTD_createDCtx();
        VERIFY(fp_.dctx);
        fp_.readPos = 0;
        fp_.output = { fp_.bufOut, 0, 0 }; // Nothing here
        fp_.input = { fp_.bufIn, 0, 0 }; // Nothing here

        // STEP 3: initialize seq
        seq_ = fastafastqzstd::kseq_init(&fp_);
        eof_ = false;
        is_open_ = true;
        ReadAhead();
    }

    /*
     * Read next SingleRead from file.
     */
    void ReadAhead() {
        VERIFY(is_open_);
        VERIFY(!eof_);
        if (fastafastqzstd::kseq_read(seq_) < 0) {
            eof_ = true;
        }
    }

    /*
     * Hidden copy constructor.
     */
    FastaFastqZstdParser(const FastaFastqZstdParser& parser);
    /*
     * Hidden assign operator.
     */
    void operator=(const FastaFastqZstdParser& parser);
};

}
