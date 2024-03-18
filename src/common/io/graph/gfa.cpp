
//***************************************************************************
//* Copyright (c) 2023-2024 SPAdes team
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

// Copyright 2022 Anton Korobeynikov

// This file is part of Bandage

// Bandage is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.

// Bandage is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.

// You should have received a copy of the GNU General Public License
// along with Bandage.  If not, see <http://www.gnu.org/licenses/>.

#include "gfa.hpp"

#include <lexy/action/parse.hpp> // lexy::parse
#include <lexy/action/trace.hpp> // lexy::trace
#include <lexy/dsl.hpp>          // lexy::dsl::*
#include <lexy/callback.hpp>     // lexy callbacks
#include <lexy/grammar.hpp>
#include <lexy/input/string_input.hpp>
#include <lexy/input/argv_input.hpp>
#include <lexy_ext/report_error.hpp>
#include <optional>

namespace gfa {

std::ostream &operator<<(std::ostream &s, const tag &t) {
    s << t.name[0] << t.name[1] << ':';
    return std::visit([&](const auto& value) -> std::ostream& { return s << value; }, t.val);
}

namespace grammar {
namespace dsl = lexy::dsl;

struct tag {
    struct tag_character : lexy::token_production {
        static constexpr auto rule = dsl::capture(dsl::ascii::alpha_digit);
        static constexpr auto value = lexy::as_string<std::string>;
    };

    struct tag_integer : lexy::token_production {
        static constexpr auto rule =
                dsl::minus_sign + dsl::integer<std::int64_t>(dsl::digits<>.no_leading_zero());
        static constexpr auto value = lexy::as_integer<std::int64_t>;
    };

    struct tag_string : lexy::token_production {
        static constexpr auto rule = dsl::identifier(dsl::ascii::print);
        static constexpr auto value = lexy::as_string<std::string>;
    };

    struct tag_float : lexy::token_production {
        static constexpr auto rule = [] {
            auto integer = dsl::if_(dsl::lit_c<'-'>) + dsl::digits<>.no_leading_zero();
            auto fraction  = dsl::lit_c<'.'> >> dsl::digits<>;
            auto exp_char = dsl::lit_c<'e'> | dsl::lit_c<'E'>;
            auto exponent = exp_char >> (dsl::lit_c<'+'> | dsl::lit_c<'-'>) + dsl::digits<>;
            return dsl::peek(dsl::lit_c<'-'> / dsl::digit<>) >>
                    dsl::position +
                    integer +
                    dsl::if_(fraction) +
                    dsl::if_(exponent) +
                    dsl::position;
        }();

        static constexpr auto value = lexy::callback<float>(
            // std::from_chars(const char*, const char*, float) is only
            // available starting from libc++ from LLVM 14 :(
            [](const char *first, const char *) { return ::atof(first); }
        );
    };

    struct tag_name : lexy::token_production {
        static constexpr auto rule = dsl::capture(dsl::token(dsl::ascii::alpha + dsl::ascii::alpha_digit));
        static constexpr auto value = lexy::as_string<std::string_view>;
    };

    static constexpr auto rule = [] {
        auto colon = dsl::lit_c<':'>;
        return dsl::p<tag_name> + colon +
        (
            dsl::capture(LEXY_LIT("A")) >> colon + dsl::p<tag_character> |
            dsl::capture(LEXY_LIT("i")) >> colon + dsl::p<tag_integer>   |
            dsl::capture(LEXY_LIT("f")) >> colon + dsl::p<tag_float>     |
            dsl::capture(LEXY_LIT("Z")) >> colon + dsl::p<tag_string>    |
            dsl::capture(LEXY_LIT("J")) >> colon + dsl::p<tag_string>    |
            dsl::capture(LEXY_LIT("H")) >> colon + dsl::p<tag_string>    |
            dsl::capture(LEXY_LIT("B")) >> colon + dsl::p<tag_string>
        );
    }();

    static constexpr auto value = lexy::callback<gfa::tag>(
        [](std::string_view name, auto type, auto val) {
            return gfa::tag{name, std::string_view{type.data(), type.size()}, val};
        });
};

auto tab = dsl::lit_c<'\t'>;

struct segment_name {
    static constexpr auto name = "segment name";
    static constexpr auto rule =
            dsl::identifier(dsl::ascii::graph - LEXY_LIT("=") - LEXY_LIT("*") - LEXY_LIT(",") - LEXY_LIT(";"),
                            dsl::ascii::graph - LEXY_LIT(",") - LEXY_LIT(";"));
    static constexpr auto value = lexy::as_string<std::string_view>;
};

struct segment_orientation {
    static constexpr auto name = "segment orientation";
    static constexpr auto rule = dsl::capture(LEXY_ASCII_ONE_OF("+-"));
    static constexpr auto value = lexy::as_string<std::string_view>;
};

struct oriented_segment {
    // Apparently we cannot use segment_name + segment_orientation as GFA grammar is context-dependent
    // Parse as full segment name and deal with possible invalid input later
    // static constexpr auto rule = dsl::capture(dsl::token(dsl::p<segment_name> + dsl::p<segment_orientation>));
    static constexpr auto rule =
            dsl::identifier(dsl::ascii::graph - LEXY_LIT("=") - LEXY_LIT("*") - LEXY_LIT(",") - LEXY_LIT(";"),
                            dsl::ascii::graph - LEXY_LIT(",") - LEXY_LIT(";"));
    static constexpr auto value = lexy::as_string<std::string_view>;
};

struct cigar_string {
    static constexpr auto name = "CIGAR string";

    static constexpr auto cigaropcode =
            LEXY_CHAR_CLASS("CIGAR opcode",
                            LEXY_LIT("M") / LEXY_LIT("I") / LEXY_LIT("D") /
                            LEXY_LIT("N") / LEXY_LIT("S") / LEXY_LIT("H") /
                            LEXY_LIT("P") / LEXY_LIT("X") / LEXY_LIT("="));

    struct cigarop : lexy::transparent_production {
        static constexpr auto name = "CIGAR operation";

        static constexpr auto rule = dsl::integer<std::uint32_t> >> dsl::capture(cigaropcode);
        static constexpr auto value = lexy::callback<gfa::cigarop>(
            [](std::uint32_t cnt, auto lexeme) {
                return gfa::cigarop{cnt, lexeme[0]};
            });
    };

    static constexpr auto rule = dsl::list(dsl::p<cigarop>);
    static constexpr auto value = lexy::as_list<std::vector<gfa::cigarop>>;
};

struct opt_tags {
    static constexpr auto rule = [] {
        auto tags = dsl::list(dsl::p<tag>, dsl::sep(tab));
        return dsl::eof | (tab >> tags + dsl::eof);
    }();
    static constexpr auto value = lexy::as_list<std::vector<gfa::tag>>;
};

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
        static constexpr auto rule = dsl::list(dsl::p<oriented_segment>, dsl::sep(dsl::comma));
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
               dsl::p<path> |
               comment |
               // Explicitly ignore all other records (though require proper tab-delimited format)
               dsl::ascii::alpha >> tab + dsl::until(dsl::newline).or_eof() |
               dsl::error<expected_gfa_record>;
     }();

    static constexpr auto value = lexy::construct<gfa::record>;
};

}; // namespace grammar


std::optional<gfa::record> parse_record(const char* line, size_t len) {
    auto result = lexy::parse<grammar::record>(lexy::string_input(line, len), lexy_ext::report_error);
    if (result.has_value())
        return std::make_optional(std::move(result.value()));

    return {};
}

};
