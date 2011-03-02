/*
 * nucl.cpp
 *
 *  Created on: 01.03.2011
 *      Author: vyahhi
 */

#include "nucl.hpp"

char complement(char c) {
	return c ^ 3;
}

char nucl(char c) {
	switch(c) {
		case 0: return 'A';
		case 1: return 'C';
		case 2: return 'G';
		case 3: return 'T';
		default: return 'N';
	}
}

char to_char(char nucl) {
	switch(nucl) {
		case 'A': return 0;
		case 'C': return 1;
		case 'G': return 2;
		case 'T': return 3;
		default: return 100;
	}
}
