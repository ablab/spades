/*
 * Simple operations and checks for nucleotide-letters
 *
 *  Created on: 01.03.2011
 *      Author: vyahhi
 */

#ifndef NUCL_HPP_
#define NUCL_HPP_

char complement(char c); // 0123 -> 3210
char nucl(char c); // 0123 -> ACGT
char denucl(char c); // ACGT -> 0123
bool is_nucl(char c); // is ACGT

#endif /* NUCL_HPP_ */
