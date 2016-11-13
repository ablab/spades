//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "single_read.hpp"
#include "delegating_reader_wrapper.hpp"

namespace io {

//todo refactor!!!
class IdSettingReaderWrapper: public DelegatingWrapper<SingleRead> {
    typedef DelegatingWrapper<SingleRead> base;
    size_t next_id_;
public:
    IdSettingReaderWrapper(base::ReadStreamPtrT reader, size_t start_id = 0) :
            base(reader), next_id_(start_id) {
    }

    /* virtual */
    IdSettingReaderWrapper& operator>>(SingleRead& read) {
        this->reader() >> read;
        read.ChangeName(ToString(next_id_++));
        return *this;
    }
};

class PrefixAddingReaderWrapper: public DelegatingWrapper<SingleRead> {
    typedef DelegatingWrapper<SingleRead> base;
    std::string prefix_;
public:
    PrefixAddingReaderWrapper(base::ReadStreamPtrT reader,
            const std::string& prefix) :
            base(reader), prefix_(prefix) {
    }

    /* virtual */
    PrefixAddingReaderWrapper& operator>>(SingleRead& read) {
        this->reader() >> read;
        read.ChangeName(prefix_ + read.name());
        return *this;
    }
};

//fixme currently leads to long stretches of ACGTACGT...
class FixingWrapper: public DelegatingWrapper<SingleRead> {
    typedef DelegatingWrapper<SingleRead> base;

    io::SingleRead MakeValid(const io::SingleRead& read) const {
        std::string str = read.GetSequenceString();
        for (size_t i = 0; i < str.length(); ++i) {
            if (!is_nucl(str[i]))
                str[i] = nucl(char(i % 4));
        }
        return io::SingleRead(read.name(), str);
    }

public:
    FixingWrapper(base::ReadStreamPtrT reader) :
            base(reader) {
    }

    /* virtual */
    FixingWrapper& operator>>(SingleRead& read) {
        this->reader() >> read;
        if (!read.IsValid()) {
            TRACE("Read " << read.name() << " was invalid. Fixing");
            read = MakeValid(read);
            VERIFY(read.IsValid());
        }
        return *this;
    }

private:
    DECL_LOGGER("FixingWrapper");
};

class NonNuclCollapsingWrapper: public DelegatingWrapper<SingleRead> {
    typedef DelegatingWrapper<SingleRead> base;

    io::SingleRead MakeValid(const io::SingleRead& read) const {
        std::string str = read.GetSequenceString();
        std::stringstream ss;
        for (size_t i = 0; i < read.size(); ++i) {
            if (is_nucl(str[i]))
                ss << str[i];
        }
        return io::SingleRead(read.name(), ss.str());
    }

public:
    NonNuclCollapsingWrapper(base::ReadStreamPtrT reader) :
            base(reader) {
    }

    /* virtual */
    NonNuclCollapsingWrapper& operator>>(SingleRead& read) {
        this->reader() >> read;
        if (!read.IsValid()) {
            TRACE("Read " << read.name() << " was invalid. Collapsing non-nucls");
            read = MakeValid(read);
            VERIFY(read.IsValid());
        }
        return *this;
    }

private:
    DECL_LOGGER("NonNuclCollapsingWrapper");
};

}
