//***************************************************************************
//* Copyright (c) 2023-2024 SPAdes team
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include <string>
#include <variant>
#include <string_view>
#include <optional>
#include <vector>
#include <algorithm>
#include <ostream>

#include <cstdint>
#include <cstdlib>
#include <cinttypes>
#include <cstdio>

namespace cigar {
    struct tag {
        char name[2];
        char type;
        std::variant<int64_t, std::string, float> val;

        template<typename T>
        tag(std::string_view n, std::string_view t, T v)
                : name{n[0], n[1]}, type(t.front()), val(std::move(v)) {}

        friend std::ostream &operator<<(std::ostream &s, const tag &t);

        void print() const {
            fprintf(stdout, "%c%c", name[0], name[1]);
            fputs(":", stdout);
            std::visit([&](const auto& value) { _print(value); }, val);
        }

      private:
        void _print(int64_t i) const {
            std::fprintf(stdout, "%c:%" PRId64, type, i);
        }

        void _print(const std::string &str) const {
            std::fprintf(stdout, "%c:%s", type, str.c_str());
        }

        void _print(float f) const {
            std::fprintf(stdout, "%c:%g", type, f);
        }
    };

    struct cigarop {
        uint32_t count : 24;
        char op : 8;

        void print() const {
            std::fprintf(stdout, "%u%c", count, op);
        }
    };

    using cigar_string = std::vector<cigarop>;

    static inline std::optional<tag>
    getTag(const char *name,
           const std::vector<tag> &tags) {
        auto res = std::find_if(tags.begin(), tags.end(),
                                [=](const tag &tag) {
                                    return (tag.name[0] == name[0] &&
                                            tag.name[1] == name[1]);
                                });
        if (res == tags.end())
            return {};

        return *res;
    }

    template<class T>
    std::optional<T> getTag(const char *name,
                            const std::vector<tag> &tags) {
        auto res = std::find_if(tags.begin(), tags.end(),
                                [=](const tag &tag) {
                                    return (tag.name[0] == name[0] &&
                                            tag.name[1] == name[1]);
                                });
        if (res == tags.end())
            return {};

        if (!std::holds_alternative<T>(res->val))
            return {};

        return std::get<T>(res->val);
    }

    std::optional<tag> parseTag(const char* line, size_t len);
}
