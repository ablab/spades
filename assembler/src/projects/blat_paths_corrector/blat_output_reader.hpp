#pragma once

#include <istream>
#include <array>
#include <vector>
#include <string>
#include <unordered_map>
#include <functional>

enum class Columns : size_t {
    match       = 0,
    mismatch    = 1,
    rep_match   = 2,
    Ns          = 3,
    Q_gap_count = 4,
    Q_gap_bases = 5,
    T_gap_count = 6,
    T_gap_bases = 7,
    strand      = 8,
    Q_name      = 9,
    Q_size      = 10,
    Q_start     = 11,
    Q_end       = 12,
    T_name      = 13,
    T_size      = 14,
    T_start     = 15,
    T_end       = 16,
    block_count = 17,
    blockSizes  = 18,
    qStarts     = 19,
    tStarts     = 20,
    TOTAL_COLUMNS_SIZE = 21
};

template<size_t N>
using Record = std::array<std::string, N>;

class TokenIterator {
    size_t from;
    size_t len;
    std::string const & str;
public:
    TokenIterator(std::string const & str)
        : from(0)
        , len(0)
        , str(str)
    {
        Step();
    }

    TokenIterator & operator ++ () {
        Step();
        return *this;
    }

    std::string operator * () const {
        return str.substr(from, len);
    }

    bool IsEnd() const noexcept {
        return len == 0;
    }

private:
    void Step() noexcept {
        from += len;
        len = 0;
        while (from < str.size() && isspace(str[from]))
            ++from;
        while (from + len < str.size() && !isspace(str[from + len]))
            ++len;
    }
};

template<size_t N>
class Records {
    using Indexes = std::array<size_t, static_cast<size_t>(Columns::TOTAL_COLUMNS_SIZE)>;
    std::vector<Record<N>> records;
    Indexes complessed_index;
public:

    class RecordProxy {
        Indexes const & complessed_index;
        Record<N>  const & record;
    public:
        RecordProxy(Indexes const & complessed_index, Record<N> const & record)
            : complessed_index(complessed_index)
            , record(record)
        {}
        
        std::string const & operator [](Columns index) const noexcept {
            auto real_index = complessed_index[static_cast<size_t>(index)];
            VERIFY(real_index != static_cast<size_t>(Columns::TOTAL_COLUMNS_SIZE));
            return record[real_index];
        }
    };

    Records(std::array<Columns, N> const & columns) {
        complessed_index.fill(static_cast<size_t>(Columns::TOTAL_COLUMNS_SIZE));
        std::array<bool, static_cast<size_t>(Columns::TOTAL_COLUMNS_SIZE)> used;
        used.fill(false);
        for (size_t i = 0; i < columns.size(); ++i) {
            if (columns[i] == Columns::TOTAL_COLUMNS_SIZE)
                throw std::string("You cannot use Columns::TOTAL_COLUMNS_SIZE in columns list");
            auto column_id = static_cast<size_t>(columns[i]);
            if (used[column_id])
                throw std::string("Column with index ") + std::to_string(column_id) + std::string(" is specified more than once");
            used[column_id] = true;
            complessed_index[column_id] = i;
        }
    }

    RecordProxy operator [] (size_t pos) const noexcept {
        return RecordProxy(complessed_index, records[pos]);
    }

    size_t size() const noexcept {
        return records.size();
    }

    void PushBack(Record<N> && record) {
        records.push_back(std::move(record));
    }

    Indexes const & GetIndexes() const noexcept {
        return complessed_index;
    }

};

template<size_t N>
using FilterType = std::function<bool(typename Records<N>::RecordProxy const &)>;

template<size_t N>
class RecordPusher {
private:
    Records<N> & records;
    FilterType<N> const & filter;
public:
    RecordPusher(Records<N> & records, FilterType<N> const & filter) 
        : records(records)
        , filter(filter)
    {}

    bool Push(std::string const & line) {
        TokenIterator it(line);
        Record<N> rd;
        size_t i = 0;
        for (size_t used = 0; i < static_cast<size_t>(Columns::TOTAL_COLUMNS_SIZE) && !it.IsEnd() && used < N; ++it, ++i) {
            if (records.GetIndexes()[i] != static_cast<size_t>(Columns::TOTAL_COLUMNS_SIZE))
                rd[records.GetIndexes()[i]] = *it;
        }
        if ((i < static_cast<size_t>(Columns::TOTAL_COLUMNS_SIZE)) != (!it.IsEnd()))
            throw std::string("blat output file is corrupted");

        auto should_be_pushed = filter(typename Records<N>::RecordProxy(records.GetIndexes(), rd));
        if (should_be_pushed)
            records.PushBack(std::move(rd));
        return should_be_pushed;
    }
};

bool GetNextNonemptyLine(std::istream & inp, std::string & result) {
    while (std::getline(inp, result)) {
        if (result.empty())
            continue;
        return true;
    }
    return false;
}

void SkipHeader(std::istream & inp) {
    std::string line;
    if (!GetNextNonemptyLine(inp, line))
        throw std::string("Empty file");
    if (line != "psLayout version 3")
        throw std::string("Sorry, unsupported blat output version");
    GetNextNonemptyLine(inp, line); // column's_names_upper_part
    GetNextNonemptyLine(inp, line); // column's_names_lower_part
    GetNextNonemptyLine(inp, line); // -------------------------
    if (line != std::string(159, '-'))
        throw std::string("blat output file is corrupted");
}

template<size_t N>
Records<N> Read(std::istream & inp,
                std::array<Columns, N> const & columns,
                FilterType<N> const & filter)
{
    Records<N> records(columns);
    RecordPusher<N> pusher(records, filter);
    std::string line;
    SkipHeader(inp);
    size_t total_lines = 0;
    size_t accepted_lines = 0;
    while (GetNextNonemptyLine(inp, line)) {
        accepted_lines += pusher.Push(line);
        ++total_lines;
    }
    INFO("Total line read: " << total_lines);
    INFO("Accepted lines: " << accepted_lines);
    return records;
}
