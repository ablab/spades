//***************************************************************************
//* Copyright (c) 2023-2024 SPAdes team
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include <lexy/action/parse.hpp> // lexy::parse
#include <lexy/input/string_input.hpp>
#include <lexy_ext/report_error.hpp>

#include "cigar.hpp"

#include "cigar.inl"

namespace cigar {
std::ostream &operator<<(std::ostream &s, const tag &t) {
    s << t.name[0] << t.name[1] << ':';
    return std::visit([&](const auto& value) -> std::ostream& { return s << value; }, t.val);
}

std::optional<tag> parseTag(const char* line, size_t len) {
    lexy::visualization_options opts;
    opts.max_lexeme_width = 35;

    auto result = lexy::parse<grammar::tag>(lexy::string_input(line, len), lexy_ext::report_error.opts(opts));
    if (result.has_value())
        return std::make_optional(result.value());

    return {};
}


}
