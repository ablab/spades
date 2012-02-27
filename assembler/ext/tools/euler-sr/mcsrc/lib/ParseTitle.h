/***************************************************************************
 * Title:          ParseTitle.h 
 * Author:         Mark Chaisson
 * Created:        2007
 * Last modified:  12/18/2008
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#ifndef PARSE_TITLE_H_
#define PARSE_TITLE_H_
#include <string>
#include <sstream>

ssize_t ParseTitle(std::string &titleIn, std::string &name);
ssize_t ExtractQuotedString(std::string &in, std::string &value);


// Parse strings of the format  "keyword=value<whitespace>"
template<typename T>
ssize_t ParseKeyword(std::string &name, std::string keyword, T& value) {
	ssize_t keywordPos;
	keywordPos = name.find(keyword);
	if (keywordPos != name.npos) {
		keywordPos += keyword.size();
		while(keywordPos < name.size() and 
					name[keywordPos] != '=')
			keywordPos++;
		// Make sure we didn't look too far.
		keywordPos++;
		if (keywordPos >= name.size() )
			return 0;

		std::stringstream namestrm(name.substr(keywordPos));
		// Find a value iwth bounda ries 
		if (!(namestrm >> value))
			return 0;
		else 
			return 1;
	}
	else {
		return 0;
	}
}
	
	
#endif
