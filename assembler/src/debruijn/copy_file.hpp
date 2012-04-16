/*
 * copy_file.hpp
 *
 *  Created on: 8 Sep 2011
 *      Author: valery
 */
#pragma once

#include <vector>
#include <boost/filesystem.hpp>

namespace fs = boost::filesystem;
namespace details
{
void copy_file(const fs::path& from_path, const fs::path& to_path);

inline void hard_link_file(const fs::path& from_path, const fs::path& to_path)
{
#if not defined(BOOST_FILESYSTEM_VERSION) or (BOOST_FILESYSTEM_VERSION == 2)
	try {
		fs::create_hard_link(from_path, to_path);
	} catch( ... ) {
		copy_file(from_path, to_path);
	}
#elif BOOST_FILESYSTEM_VERSION == 3
	try {
		boost::filesystem3::create_hard_link(from_path, to_path);
	} catch( ... ) {
		copy_file(from_path, to_path);
	}
#else
	BOOST_STATIC_ASSERT(false && "BOOST_FILESYSTEM_VERSION defined, but has value different from 2 or 3");
#endif
}

inline void copy_file(const fs::path& from_path, const fs::path& to_path)
{
#if not defined(BOOST_FILESYSTEM_VERSION) or (BOOST_FILESYSTEM_VERSION == 2)
	fs::copy_file(from_path, to_path);
#elif BOOST_FILESYSTEM_VERSION == 3
	boost::filesystem3::copy(from_path, to_path);
#else
	BOOST_STATIC_ASSERT(false && "BOOST_FILESYSTEM_VERSION defined, but has value different from 2 or 3");
#endif
}
}


typedef std::vector<fs::path> files_t;

files_t files_by_prefix(fs::path const& p)
{
	//using namespace details;
	files_t files;

	fs::path folder(p.parent_path());
	std::string prefix = p.filename().c_str();

	for (auto it = fs::directory_iterator(folder), end = fs::directory_iterator(); it != end; ++it)
	{
		if (is_regular_file(*it) && boost::starts_with(it->path().filename().c_str(), prefix))
			files.push_back(*it);
	}

	return files;
}

void copy_files_by_prefix(files_t const& files, fs::path const& to_folder)
{
	for (size_t i = 0; i != files.size(); ++i)
	{
		fs::path f(files[i]);
		files_t  files_to_copy = files_by_prefix(f);

		for (auto it = files_to_copy.begin(); it != files_to_copy.end(); ++it)
			details::copy_file(*it, to_folder / it->filename());
	}
}

void link_files_by_prefix(files_t const& files, fs::path const& to_folder)
{
	for (size_t i = 0; i != files.size(); ++i)
	{
		fs::path f(files[i]);
		files_t  files_to_copy = files_by_prefix(f);

		for (auto it = files_to_copy.begin(); it != files_to_copy.end(); ++it)
			details::hard_link_file(*it, to_folder / it->filename());
	}
}

void copy_files_by_ext(fs::path const& from_folder, fs::path const& to_folder, std::string const& ext, bool recursive)
{
	for (auto it = fs::directory_iterator(from_folder), end = fs::directory_iterator(); it != end; ++it)
    {
	    if (recursive && is_directory(it->path()))
        {
		    fs::path subdir = to_folder / it->path().filename();
		    make_dir(subdir);
		    copy_files_by_ext(it->path(), subdir, ext, recursive);
        }

	    if (it->path().extension() == ext)
	    	details::copy_file(*it, to_folder / it->path().filename());
    }
}
