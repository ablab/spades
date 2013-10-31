#pragma once
#include "standard_base.hpp"
#include "io/delegating_reader_wrapper.hpp"

namespace cap {

class RepeatCroppingReader : public io::DelegatingReaderWrapper<io::SingleRead> {
    typedef io::DelegatingReaderWrapper<io::SingleRead> base;
    vector<pair<size_t, size_t>> coordinates_ladder_;
public:
    RepeatCroppingReader(io::IReader<io::SingleRead>& reader) : base(reader) {

    }

    RepeatCroppingReader& operator>>(io::SingleRead& read) {
        base::operator >>(read);
        coordinates_ladder_.clear();
        string orig_string = read.GetSequenceString();
        string orig_qual = read.GetQualityString();
        string cropped = "";
        string cropped_qual = "";
        coordinates_ladder_.push_back(make_pair(0, 0));
        for (size_t coord = 0; coord < orig_string.size(); ++coord) {
            if (isupper(orig_string[coord])) {
                cropped += orig_string[coord];
                cropped_qual += orig_qual[coord];
            }
            if (coord > 0 && coord < orig_string.size() - 1
                    && (isupper(orig_string[coord]) ^ isupper(orig_string[coord + 1]))) {
                coordinates_ladder_.push_back(make_pair(coord, cropped.size() - 1));
            }
        }
        coordinates_ladder_.push_back(make_pair(orig_string.size(), cropped.size()));
        read = io::SingleRead(read.name(), cropped, cropped_qual);
        return *this;
    }

    vector<pair<size_t, size_t>> coordinates_ladder() const {
        return coordinates_ladder_;
    }
};

}
