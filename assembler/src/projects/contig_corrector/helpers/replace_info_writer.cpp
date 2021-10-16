#include "replace_info_writer.hpp"
#include <iomanip>
#include <cassert>

namespace helpers {

namespace {

const char* ToString(RangeType type) {
    #define CASE(w) case RangeType::w: return #w; break
    switch (type) {
        CASE(origin);
        CASE(edge);
        CASE(path);
        default: assert(false && "unreachable");
    }
    #undef CASE
}

template<class T>
void Print(std::ofstream & out, const char * name, const char * type, T from, T len) {
    using std::setw;
    out << setw(15) << name << ' ' << setw(10) << type << ' ' << setw(10) << from << ' ' << setw(10) << len << '\n';
}

} // namespace


ReplaceInfoWriter::ReplaceInfoWriter(std::string file)
    : out(file)
{
    if (!out.is_open())
        throw "Cannot create '" + file + "' file";
    Print(out, "name", "seq_type", "from", "len");
}

void ReplaceInfoWriter::Write(std::string const & name, RangeType type, unsigned long long from, unsigned long long len) {
    Print(out, name.c_str(), ToString(type), from, len);
}

std::unique_ptr<ReplaceInfoWriter> ReplaceInfoWriter::stream;

void ReplaceInfoWriter::SetStream(std::string file) {
    stream.reset(new ReplaceInfoWriter(std::move(file)));
}

} // namespace helpers
