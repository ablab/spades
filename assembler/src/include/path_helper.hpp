#pragma once

#include <boost/tokenizer.hpp>

#include <sys/types.h>
#include <sys/stat.h>
#include <dirent.h>
#include <unistd.h>

namespace path
{

inline bool make_dir(string const& folder)
{
    return mkdir(folder.c_str(), 0755) == 0;
}

inline void remove_dir(string const& folder)
{
    DIR *dp;
    if ((dp = opendir(folder.c_str())) == NULL)
        throw std::runtime_error("can not open folder " + folder);

    struct dirent *dirp;
    while ((dirp = readdir(dp)) != NULL)
    {
        string full_path = folder + "/" + dirp->d_name;

        if (dirp->d_type == DT_DIR)
        {
            if (string("." ) != dirp->d_name &&
                string("..") != dirp->d_name)
            {
                remove_dir(full_path);
            }
        }
        else
            remove(full_path.c_str());
    }

    closedir(dp);
    remove  (folder.c_str());
}

inline bool is_regular_file(string const& path)
{
    struct stat st;
    return (stat(path.c_str(), &st) == 0) && (S_ISREG(st.st_mode));
}

inline string append_path(string const& prefix, string const& suffix)
{
    string delimiter = "";

    if (!boost::ends_with  (prefix, "/") &&
        !boost::starts_with(suffix, "/") &&
        !prefix.empty())
    {
        delimiter = "/";
    }

    return prefix + delimiter + suffix;
}

inline string current_dir()
{
    char* cwd = get_current_dir_name();
    string result = cwd;

    free(cwd);
    return result;
}

inline void make_full_path(string& path)
{
    if (!boost::starts_with(path, "/")) // relative path
        path = append_path(current_dir(), path);
}

inline string filename(string const& path)
{
    size_t pos = path.find_last_of('/');
    return pos != string::npos ? path.substr(pos + 1) : path;
}

inline string basename(string const& path)
{
    size_t slash       = path.find_last_of('/');
    size_t after_slash = slash == string::npos ? 0 : slash + 1;

    size_t dot   = path.find_last_of('.');
    if (dot < after_slash)
        dot = string::npos;

    return path.substr(after_slash, dot - after_slash);
}

inline string extension(string const& path)
{
    size_t slash       = path.find_last_of('/');
    size_t after_slash = slash == string::npos ? 0 : slash + 1;
    size_t dot	 	   = path.find_last_of('.');

    if (dot < after_slash   ||
        dot == string::npos ||
        dot + 1 == path.size())

        return string();

    return path.substr(dot);
}

inline string parent_path(string path)
{
    make_full_path(path);
    size_t slash_pos = path.find_last_of('/');

    return slash_pos == 0
        ? string("/")
        : path.substr(0, slash_pos);
}

// doesn't support symlinks
inline string resolve(string const& path)
{
    typedef boost::char_delimiters_separator<char> 	separator_t;
    typedef boost::tokenizer<separator_t>			tokenizer_t;

    tokenizer_t tok(path, separator_t(false, "", "/"));

    string result = "/";
    for (auto it = tok.begin(); it != tok.end(); ++it)
    {
        if(*it == "..")
            result = parent_path(result);

        else if(*it == ".")
            ;// Ignore

        else
            // Just cat other path entries
            result = append_path(result, *it);
    }

    return result;
}

inline string make_relative_path(string p, string base = current_dir())
{
    p    = resolve(p);
    base = resolve(base);

    string pp = parent_path(p);

    typedef boost::char_delimiters_separator<char> 	separator_t;
    typedef boost::tokenizer<separator_t>			tokenizer_t;

    tokenizer_t pp_tok  (pp  , separator_t(false, "", "/"));
    tokenizer_t base_tok(base, separator_t(false, "", "/"));

    auto i = pp_tok  .begin();
    auto j = base_tok.begin();

    while (i != pp_tok.end() && j != base_tok.end() && *i == *j)
    {
        ++i;
        ++j;
    }

    string result;
    for ( ; j != base_tok.end(); ++j)
        result = append_path("..", result);

    for ( ; i != pp_tok.end(); ++i)
        result = append_path(result, *i);

    return append_path(result, filename(p));
}

typedef vector<string> files_t;

}
