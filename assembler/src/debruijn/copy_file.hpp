/*
 * copy_file.hpp
 *
 *  Created on: 8 Sep 2011
 *      Author: valery
 */
#pragma once

typedef std::vector<fs::path> files_t;

files_t files_by_prefix(fs::path const& p)
{
	files_t files;

	fs::path folder(p.parent_path());
	std::string prefix = p.filename();

	for (auto it = fs::directory_iterator(folder), end = fs::directory_iterator(); it != end; ++it)
	{
		if (is_regular_file(*it) && boost::starts_with(it->filename(), prefix))
			files.push_back(*it);
	}

	return files;
}


void copy_files(files_t const& files)
{
	for (size_t i = 0; i != files.size(); ++i)
	{
		fs::path f(files[i]);
		files_t  files_to_copy = files_by_prefix(f);
		fs::path to_folder = cfg::get().output_saves;

		for (auto it = files_to_copy.begin(); it != files_to_copy.end(); ++it)
			copy_file(*it, to_folder / it->filename());
	}
}
