/***************************************************************************
 * Title:          ParseTitle.cpp 
 * Author:         Mark Chaisson
 * Created:        2007
 * Last modified:  12/18/2008
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#include "ParseTitle.h"


ssize_t ParseTitle(std::string &titleIn, std::string &name) {
	name = titleIn;
	if (name.size() == 0) return 0;

	if (name.c_str()[0] == '>')
		name.replace(0,1, "");
	ssize_t len = name.size();
	if (name[len-1] == '\r')
		name.erase(len-1,len);

	ssize_t blankIndex = name.find(" ");
	if (blankIndex != name.npos) {
		name = name.substr(0, blankIndex);
	}


	return 1;
}


ssize_t ExtractQuotedString(std::string &in, std::string &value) {
	ssize_t i;
	ssize_t begin = -1;
	ssize_t end = -1;
	for (i =0; i < in.size(); i++) {
		if (in[i] == '"') {
			begin = i; break;
		}
	}
	for (i++; i < in.size(); i++ ) {
		if (in[i] == '"' and in[i-1] != '\\') {
			end = i+1;
			break;
		}
	}
	value = "";

	if (begin >= 0) 
		value = in.substr(begin+1, end-begin-2);
	return 1;
}

