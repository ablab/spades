//***************************************************************************
//* Copyright (c) 2023-2024 SPAdes team
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include <lexy/dsl.hpp>          // lexy::dsl::*
#include <lexy/callback.hpp>     // lexy callbacks
#include <lexy/grammar.hpp>

#include <string>

namespace cigar::grammar {
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
                auto integer = dsl::if_(dsl::lit_c < '-' > ) + dsl::digits<>.no_leading_zero();
                auto fraction = dsl::lit_c < '.' > >> dsl::digits<>;
                auto exp_char = dsl::lit_c < 'e' > | dsl::lit_c<'E'>;
                auto exponent = exp_char >> (dsl::lit_c < '+' > | dsl::lit_c < '-' > ) + dsl::digits<>;
                return dsl::peek(dsl::lit_c < '-' > / dsl::digit<>) >>
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
            static constexpr auto name = "tag name";

            static constexpr auto rule = dsl::capture(dsl::token(dsl::ascii::alpha + dsl::ascii::alpha_digit));
            static constexpr auto value = lexy::as_string<std::string_view>;
        };

        struct invalid_tag_type {
            static constexpr auto name = "invalid tag type";
        };

        static constexpr auto rule = [] {
            auto colon = dsl::lit_c<':'>;
            return dsl::p<tag_name> >> colon +
                   (
                           dsl::capture(LEXY_LIT("A")) >> colon + dsl::p < tag_character > |
                           dsl::capture(LEXY_LIT("i")) >> colon + dsl::p < tag_integer > |
                           dsl::capture(LEXY_LIT("f")) >> colon + dsl::p < tag_float > |
                           dsl::capture(LEXY_LIT("Z")) >> colon + dsl::p < tag_string > |
                           dsl::capture(LEXY_LIT("J")) >> colon + dsl::p < tag_string > |
                           dsl::capture(LEXY_LIT("H")) >> colon + dsl::p < tag_string > |
                           dsl::capture(LEXY_LIT("B")) >> colon + dsl::p < tag_string > |
                           dsl::error<invalid_tag_type>
                   );
        }();

        static constexpr auto value = lexy::callback<cigar::tag>(
                [](std::string_view name, auto type, auto val) {
                    return cigar::tag{name, std::string_view{type.data(), type.size()}, val};
                });
    };

    struct cigar_string {
        static constexpr auto name = "CIGAR string";

        static constexpr auto cigaropcode =
                LEXY_CHAR_CLASS("CIGAR opcode",
                                LEXY_LIT("M") / LEXY_LIT("I") / LEXY_LIT("D") /
                                LEXY_LIT("N") / LEXY_LIT("S") / LEXY_LIT("H") /
                                LEXY_LIT("P") / LEXY_LIT("X") / LEXY_LIT("=")) / LEXY_LIT("J");

        struct cigarop : lexy::transparent_production {
            static constexpr auto name = "CIGAR operation";

            static constexpr auto rule =
                    dsl::period |
                    dsl::integer<std::uint32_t> >> dsl::capture(cigaropcode);
            static constexpr auto value = lexy::callback<cigar::cigarop>(
                    []() { return cigar::cigarop{0, 0}; },
                    [](std::uint32_t cnt, auto lexeme) {
                        return cigar::cigarop{cnt, lexeme[0]};
                    });
        };

        static constexpr auto rule = dsl::list(dsl::p<cigarop>);
        static constexpr auto value = lexy::as_list<std::vector<cigar::cigarop>>;
    };

    static constexpr auto tab = dsl::lit_c<'\t'>;

    struct opt_tags {
        static constexpr auto name = "tags";

        static constexpr auto rule = [] {
            auto tags = dsl::list(dsl::p<tag>, dsl::trailing_sep(tab));
            return dsl::eof | (tab >> tags + dsl::eof);
        }();
        static constexpr auto value = lexy::as_list<std::vector<cigar::tag>>;
    };
}
