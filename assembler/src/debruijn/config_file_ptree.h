/*
 * config_file_ptree.h
 *
 *  Created on: Aug 9, 2011
 *      Author: Alexey.Gurevich
 */

#ifndef CONFIG_FILE_PTREE_H_
#define CONFIG_FILE_PTREE_H_

#include <boost/property_tree/ptree.hpp>
#include <string>
#include <fstream>

using boost::property_tree::ptree;
using boost::property_tree::ptree_bad_path;
using std::string;

void test_cfg_ptree();

// Meyers' singleton
class config
{
public:
	static config& get_instance(string const& filename,
            					string const& delimiter = "=",
            					string const& comment = "#")
	{
		static config* cnf = new config(filename, delimiter, comment);
		return *cnf;
	}

	template <class T>
	T get(string const& key)
	{
		try
		{
			return config_ptree_.get<T>(key);
		}
		catch (ptree_bad_path const& pbp)
		{
			throw key_not_found(key);
		}
	}

private:
	config(string const& filename,
           string const& delimiter,
           string const& comment)
	{
		std::ifstream in(filename.c_str());
		if (!in) throw file_not_found(filename);

		while (in)
		{
			// Read an entire line at a time
			string line;
			std::getline(in, line);

			// Ignore comments
			line = line.substr(0, line.find(comment));
			if (line.length() == 0)
				continue;

			// Parse the line if it contains a delimiter
			size_t delimPos = line.find(delimiter);
			if (delimPos != string::npos) // if found
			{
				// Extract the key
				string key = line.substr(0, delimPos);
				string value = line.substr(delimPos + delimiter.length());

				// Store key and value
				trim(key);
				trim(value);
				config_ptree_.put(key, value);
			}
		}
	}

	void trim(string& s)
	{
		// Remove leading and trailing whitespace
		static const char whitespace[] = " \n\t\v\r\f";
		s.erase( 0, s.find_first_not_of(whitespace) );
		s.erase( s.find_last_not_of(whitespace) + 1U );
	}

	ptree config_ptree_;

// Exception types
public:
	struct file_not_found
	{
		string filename_;
		file_not_found(string const& filename = string() )
			: filename_(filename) {}
	};

	struct key_not_found
	{
		string key_;
		key_not_found(string const& key = string() )
			: key_(key) {}
	};
};

#endif /* CONFIG_FILE_PTREE_H_ */
