//***************************************************************************
//* Copyright (c) 2023-2024 SPAdes team
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include "gfa.hpp"
#include "cigar.hpp"

#include <lexy/action/parse.hpp> // lexy::parse
#include <lexy/action/trace.hpp> // lexy::trace
#include <lexy/dsl.hpp>          // lexy::dsl::*
#include <lexy/callback.hpp>     // lexy callbacks
#include <lexy/grammar.hpp>
#include <lexy/input/string_input.hpp>
#include <lexy/input/argv_input.hpp>
#include <lexy_ext/report_error.hpp>
#include <optional>

#include "cigar.inl"

namespace gfa {

namespace grammar {
namespace dsl = lexy::dsl;

using cigar::grammar::tab;

struct segment_name {
    static constexpr auto name = "segment name";
    static constexpr auto rule =
            dsl::identifier(dsl::ascii::graph - LEXY_LIT("=") - LEXY_LIT("*") - LEXY_LIT(",") - LEXY_LIT(";"),
                            dsl::ascii::graph - LEXY_LIT(",") - LEXY_LIT(";") - LEXY_LIT("<") - LEXY_LIT(">"));
    static constexpr auto value = lexy::as_string<std::string_view>;
};

struct segment_orientation {
    static constexpr auto name = "segment orientation";
    static constexpr auto rule = dsl::capture(LEXY_ASCII_ONE_OF("+-"));
    static constexpr auto value = lexy::as_string<std::string_view>;
};

struct segment_distance {
    static constexpr auto rule =
            dsl::sign + dsl::integer<std::int64_t>(dsl::digits<>.no_leading_zero());
    static constexpr auto value = lexy::as_integer<std::int64_t>;
};

struct oriented_segment {
    // Apparently we cannot use segment_name + segment_orientation as GFA grammar is context-dependent
    // Parse as full segment name and deal with possible invalid input later
    // static constexpr auto rule = dsl::capture(dsl::token(dsl::p<segment_name> + dsl::p<segment_orientation>));
    static constexpr auto rule =
            dsl::identifier(dsl::ascii::graph - LEXY_LIT("=") - LEXY_LIT("*") - LEXY_LIT(",") - LEXY_LIT(";"),
                            dsl::ascii::graph - LEXY_LIT(",") - LEXY_LIT(";") - LEXY_LIT("<") - LEXY_LIT(">"));
    static constexpr auto value = lexy::as_string<std::string_view>;
};

using cigar::grammar::cigar_string;
using cigar::grammar::opt_tags;

// Header
// ======
// Required fields:
// Column    Field        Type        Regexp    Description
// 1         RecordType   Character   H         Record type
struct header {
    static constexpr auto rule = LEXY_LIT("H") >> dsl::p<opt_tags>;
    static constexpr auto value = lexy::construct<gfa::header>;
};

// Segment
// =======
// Required fields:
// Column    Field       Type        Regexp               Description
// 1         RecordType  Character   S                    Record type
// 2         Name        String      [!-)+-<>-~][!-~]*    Segment name
// 3         Sequence    String      \*|[A-Za-z=.]+       Optional nucleotide sequence
struct segment {
    static constexpr auto name = "GFA segment";

    static constexpr auto rule =
            LEXY_LIT("S") >>
            tab + dsl::p<segment_name> +
            tab + (LEXY_LIT("*") |
                   dsl::identifier(dsl::ascii::alpha)) +
            dsl::p<opt_tags>;
    static constexpr auto value = lexy::construct<gfa::segment>;
};

// Link line
// =========
// Required fields:
// Column   Field        Type        Regexp                   Description
// 1        RecordType   Character   L                        Record type
// 2        From         String      [!-)+-<>-~][!-~]*        Name of segment
// 3        FromOrient   String      +|-                      Orientation of From segment
// 4        To           String      [!-)+-<>-~][!-~]*        Name of segment
// 5        ToOrient     String      +|-                      Orientation of To segment
// 6        Overlap      String      \*|([0-9]+[MIDNSHPX=])+  Optional CIGAR string describing overlap
// The Overlap field is optional and can be *, meaning that the CIGAR string is not specified.
struct link {
    static constexpr auto name = "GFA link line";

    static constexpr auto rule =
            LEXY_LIT("L") >>
            tab + dsl::p<segment_name> +
            tab + dsl::p<segment_orientation> +
            tab + dsl::p<segment_name> +
            tab + dsl::p<segment_orientation> +
            tab + (LEXY_LIT("*") | dsl::p<cigar_string>) +
            dsl::p<opt_tags>;
    static constexpr auto value = lexy::construct<gfa::link>;
};

// Jump line
// =============
// Required fields:
// Column   Field        Type        Regexp                   Description
// 1        RecordType   Character   J                        Record type
// 2        From         String      [!-)+-<>-~][!-~]*        Name of segment
// 3        FromOrient   String      +|-                      Orientation of From segment
// 4        To           String      [!-)+-<>-~][!-~]*        Name of segment
// 5        ToOrient     String      +|-                      Orientation of To segment
// 6        Distance     String      \*|[-+]?[0-9]+           Optional estimated distance between the segments
struct gaplink {
    static constexpr auto name = "GFA jump line";

    static constexpr auto rule =
            LEXY_LIT("J") >>
            tab + dsl::p<segment_name> +
            tab + dsl::p<segment_orientation> +
            tab + dsl::p<segment_name> +
            tab + dsl::p<segment_orientation> +
            tab + (LEXY_LIT("*") | (dsl::else_ >> dsl::p<segment_distance>)) +
            dsl::p<opt_tags>;
    static constexpr auto value = lexy::construct<gfa::gaplink>;
};

// Path line
// =========
// Required fields
// Column   Field            Type        Regexp                   Description
// 1        RecordType       Character   P                        Record type
// 2        PathName         String      [!-)+-<>-~][!-~]*        Path name
// 3        SegmentNames     String      [!-)+-<>-~][!-~]*        A comma-separated list of segment names and orientations
// 4        Overlaps         String      \*|([0-9]+[MIDNSHPX=])+  Optional comma-separated list of CIGAR strings
struct path {
    static constexpr auto name = "GFA path record";

    struct segments {
        static constexpr auto rule = dsl::list(dsl::p<oriented_segment>, dsl::sep(dsl::comma | dsl::semicolon));
        static constexpr auto value = lexy::as_list<std::vector<std::string_view>>;
    };

    struct overlaps {
        static constexpr auto rule = dsl::list(dsl::p<cigar_string>, dsl::sep(dsl::comma));
        static constexpr auto value = lexy::as_list<std::vector<std::vector<gfa::cigarop>>>;
    };

    static constexpr auto rule =
            LEXY_LIT("P") >>
            tab + dsl::p<segment_name> +
            tab + dsl::p<segments> +
            tab + (LEXY_LIT("*") | dsl::p<overlaps>) +
            dsl::p<opt_tags>;
    static constexpr auto value = lexy::construct<gfa::path>;
};

// Walk line
// =========
// Required fields
// Column  Field        Type        Regexp              Description
// 1       RecordType   Character   W                   Record type
// 2       SampleId     String      [!-)+-<>-~][!-~]*   Sample identifier
// 3       HapIndex     Integer     [0-9]+              Haplotype index
// 4       SeqId        String      [!-)+-<>-~][!-~]*   Sequence identifier
// 5       SeqStart     Integer     \*\|[0-9]+          Optional Start position
// 6       SeqEnd       Integer     \*\|[0-9]+          Optional End position (BED-like half-close-half-open)
// 7       Walk         String      ([><][!-;=?-~]+)+   Walk
struct walk {
    static constexpr auto name = "GFA walk record";

    struct sample_id : public segment_name {
        static constexpr auto name = "sample id";
    };

    struct sequence_id : public segment_name {
        static constexpr auto name = "sequence id";
    };


    struct oriented_wsegment {
        static constexpr auto rule = dsl::capture(dsl::token(LEXY_ASCII_ONE_OF("<>") + dsl::p<segment_name>));
        static constexpr auto value = lexy::as_string<std::string_view>;
    };

    struct wsegments {
        static constexpr auto rule = dsl::list(dsl::p<oriented_segment>);
        static constexpr auto value = lexy::as_list<std::vector<std::string_view>>;
    };

    struct opt_uint64_t {
        static constexpr auto rule = LEXY_LIT("*") | dsl::integer<std::uint64_t>;
        static constexpr auto value =
                lexy::callback<gfa::walk::opt_uint64_t>(
                    [](uint64_t val) { return val; },
                    []() { return std::optional<std::uint64_t>(); });
    };

    static constexpr auto rule =
            LEXY_LIT("W") >>
            tab + dsl::p<sample_id> + // SampleId
            tab + dsl::integer<unsigned> + // HapIndex
            tab + dsl::p<sequence_id> + // SeqId
            tab + dsl::p<opt_uint64_t> + // SeqStart
            tab + dsl::p<opt_uint64_t> + // SeqEnd
            tab + dsl::p<wsegments> + // Walk
            dsl::p<opt_tags>;
    static constexpr auto value = lexy::construct<gfa::walk>;
};


struct record {
    static constexpr auto name = "GFA record";

    struct expected_gfa_record {
        static LEXY_CONSTEVAL auto name() {
            return "expected GFA record type";
        }
    };

    static constexpr auto rule = [] {
        auto comment = dsl::lit_c<'#'> >> dsl::until(dsl::newline).or_eof();

        return dsl::p<header> |
               dsl::p<segment> |
               dsl::p<link> |
               dsl::p<gaplink> |
               dsl::p<path> |
               dsl::p<walk> |
               comment |
               // Explicitly ignore all other records (though require proper tab-delimited format)
               dsl::ascii::alpha >> tab + dsl::until(dsl::newline).or_eof() |
               dsl::error<expected_gfa_record>;
     }();

    static constexpr auto value = lexy::construct<gfa::record>;
};

}; // namespace grammar


std::optional<gfa::record> parseRecord(const char* line, size_t len) {
    lexy::visualization_options opts;
    opts.max_lexeme_width = 35;

    auto result = lexy::parse<grammar::record>(lexy::string_input(line, len), lexy_ext::report_error.opts(opts));
    if (result.has_value())
        return std::make_optional(result.value());

    return {};
}

};
