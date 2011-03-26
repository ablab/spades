/*
 * bayes_quality.cpp
 *
 *  Created on: 23.03.2011
 *      Author: snikolenko
 */

#include "bayes_quality.hpp"
#include <algorithm>
#include <numeric>
#include <iostream>
#include <iomanip>
#include <math.h>

using namespace std;

LOGGER("b");

namespace bayes_quality {

double BayesQualityGenome::ReadBQInt(QRead r) {

	// compute total quality of the read
	totalq_ = accumulate(r.second.begin(), r.second.end(), 0);

	// initialize temporary vectors of possible start positions
	qvsize_ = genome_.size() - r.first.size() + 1 + INS;
	for (size_t i=0; i < INS+DEL+1; ++i) {
		fill(qv[i]->begin(), qv[i]->end(), 0);
	}
	
	INFO(r.first.str());
	// INFO(genome_.str());
	
	// dynamic programming
	for (size_t i = 0; i < r.first.size()+1; ++i) {
		for (size_t j = 0; j < qvsize_; ++j) {
			if ( i + j > genome_.size() ) { // genome is over
				for (size_t s=0; s < r.first.size()-i+1; ++s) qv[s]->at(j) = VERYLARGEQ;
				continue;
			}
			else if ( i == 0 ) qv[0]->at(j) = 0;
			else if ( r.first[i-1] != genome_[i+j-1] ) qv[0]->at(j) += r.second[i-1];
			
			// processing inserts
			for (size_t s=1; s <= INS; ++s) {
				if (i == 0) { qv[s]->at(j) = INSERT_Q*s; continue; } // first column
				if (i+s > r.first.size()) { continue; } // the read is over, too many inserts
				
				// case 1: quality value accumulated so far with s inserts + current mismatch value (down-right)
				int curq = qv[s]->at(j) + (r.first[i-1+s] != genome_[i-1+j]) * r.second[i-1+s];
				
				// case 2: quality value accumulated for s-1 inserts + insert value (down)
				int prevq = qv[s-1]->at(j) + INSERT_Q;
				
				//if ((s==1 || s==2) && j==13) cout << "i=" << i << "  " << j << " " << s << " "  << nucl(r.first[i-1+s]) << " " << curq << " " << prevq << " " << " " << " " << nucl(genome_[i-1+j]) << "\n";

				// to get the result, we add up the corresponding values and take the logarithm
				qv[s]->at(j) = AddTwoQualityValues(curq, prevq);
			}
			
			// processing deletes
			for (size_t s=INS+1; s < INS+DEL+1; ++s) {
				if (i == 0) { qv[s]->at(j) = DELETE_Q*(s-INS); continue; } // first column
				if (i+j+(s-INS) > genome_.size()) { qv[s]->at(j) = VERYLARGEQ; continue; } // the genome is over, too many deletes

				// case 1: quality value accumulated so far with s-INS deletes + current mismatch value (down-right)
				int curq = qv[s]->at(j) + (r.first[i-1] != genome_[i-1+j+s-INS]) * r.second[i-1];

				// case 2: quality value accumulated for s-INS-1 deletes + delete value (right)
				int prevq = (s==INS+1) ? (qv[0]->at(j) +   DELETE_Q) 
									   : (qv[s-1]->at(j) + DELETE_Q);
									   
				// to get the result, we add up the corresponding values and take the logarithm
				qv[s]->at(j) = AddTwoQualityValues(curq, prevq);
			}
		}
		// debug output
		/*cout << "---  i=" << i << "  " << (i < r.first.size() ? nucl(r.first[i]) : '-') << " ----\n";
		for (size_t s=0; s < INS+DEL+1; ++s) {
			ostringstream os; os.setf(ios::fixed, ios::floatfield); os.precision(1);
			for (size_t p = 0; p < qvsize_; ++p) os << setw(6) << qv[s]->at(p) << " ";
			os << "\n";
			cout << os.str();
		}*/
	}
		
	bestq_ = totalq_;
	inserts_ = 0; deletes_ = 0;
	index_ = 0;

	// now add them all up
	double res = 0;
	for (size_t i=0; i < INS+DEL+1; ++i) {
		for (size_t j = 0; j < qvsize_; ++j) {
			res += pow(10, -(qv[i]->at(j))/10.0);
			if (qv[i]->at(j) < bestq_) {
				bestq_ = (int)(qv[i]->at(j));
				index_ = j;
				if (i <= INS) { inserts_ = i; deletes_ = 0; }
				else { inserts_ = 0; deletes_ = i-INS; }
			}
		}
	}
	
	// find and pretty print the best match
	// to do so, we run dynamic programming once again (since we haven't been storing the entire matrix)
	vector<QVector *> m(r.first.size()+1);
	for (size_t i = 0; i < r.first.size()+1; ++i) {
		m[i] = new QVector(r.first.size() + deletes_ + 1);
		fill(m[i]->begin(), m[i]->end(), 0);
	}
	// first row
	for (size_t j = 0; j < r.first.size() + deletes_ + 1; ++j) {
		if (j > deletes_) { m[0]->at(j) = VERYLARGEQ; }
		else {
			m[0]->at(j) = DELETE_Q*j;
		}
	}
	// first column
	for (size_t i = 1; i < r.first.size() + 1; ++i) {
		if (i > inserts_) { m[i]->at(0) = VERYLARGEQ; }
		else {
			m[i]->at(0) = INSERT_Q*i;
		}
	}
	// other elements
	for (size_t i = 1; i < r.first.size() + 1; ++i) {
		for (size_t j = 1; j < r.first.size() + deletes_ + 1; ++j) {
			//cout << i << " " << j << " " << index_+j;
			if (i > inserts_ + j) { m[i]->at(j) = VERYLARGEQ; continue; }
			if (j > deletes_ + i) { m[i]->at(j) = VERYLARGEQ; continue; }
			m[i]->at(j) = AddThreeQualityValues(
				(index_+j <= genome_.size()) ? (m[i-1]->at(j-1) + (genome_[index_+j-1] != r.first[i-1])*r.second[i-1]) : VERYLARGEQ,
				              (deletes_ > 0) ? (m[i]->at(j-1)   + DELETE_Q)                                            : VERYLARGEQ,
				              (inserts_ > 0) ? (m[i-1]->at(j)   + INSERT_Q)                                            : VERYLARGEQ );
			//cout << " mat\n";
		}
	}
	// debug output
	/*for (size_t i = 0; i < r.first.size()+1; ++i) {
		ostringstream os; os.setf(ios::fixed, ios::floatfield); os.precision(0); os.width(4); os.fill(' ');
		for (size_t j = 0; j < r.first.size() + deletes_+1; ++j) os << setw(4) << m[i]->at(j) << " ";
		os << "\n";
		cout << os.str();
	}*/
	match_.clear();
	ostringstream osm, osmr, osmp;
	size_t ifirst_ = r.first.size();
	size_t j = min_element(m[ifirst_]->begin(), m[ifirst_]->end()) - m[ifirst_]->begin();
	for (size_t i = ifirst_; i >= 1; --i) {
		// go up
		if ((inserts_ > 0) && (j == 0 || (m[i-1]->at(j) < m[i-1]->at(j-1)) && (m[i-1]->at(j) < m[i]->at(j-1)) )) {
			match_.push_back( 2 ); osm << '-'; osmp << ' '; osmr << nucl(r.first[i-1]); continue;
		}
		// or left
		if ((deletes_ > 0) && ((m[i]->at(j-1) < m[i-1]->at(j) && m[i]->at(j-1) < m[i-1]->at(j-1)) )) {
			match_.push_back( 3 ); osm << nucl(genome_[j-1+index_]); osmp << ' '; osmr << '-'; ++i; --j; continue;
		}
		// or up-left
		match_.push_back( (int)(genome_[j-1+index_] == r.first[i-1]) );
		osm << nucl(genome_[j-1+index_]); osmr << nucl(r.first[i-1]);
		if (genome_[j-1+index_] == r.first[i-1]) osmp << "|"; else osmp << " "; --j;
	}
	// the first character is always here
	
	reverse(match_.begin(), match_.end());
	matchreadstr_ = osmr.str(); reverse(matchreadstr_.begin(), matchreadstr_.end());
	matchstr_ = osm.str(); reverse(matchstr_.begin(), matchstr_.end());
	matchprettystr_ = osmp.str(); reverse(matchprettystr_.begin(), matchprettystr_.end());
	
	for (size_t i = 0; i < r.first.size(); ++i) {
		m[i]->clear(); delete m[i];
	}
	// return the resulting likelihood
	return res;
}


float BayesQualityGenome::AddTwoQualityValues(float curq, float prevq) {
	// TODO: rewrite correctly
	return min(curq, prevq);
	
	// there are a few shortcuts here
	/*if ( curq < prevq - 1 ) return curq;
	else if ( curq > prevq + 1 ) return prevq;
	else if ( curq == prevq + 1 ) return prevq - TENLOGELEVENTENTH;
	else if ( curq == prevq - 1 ) return curq - TENLOGELEVENTENTH;
	else if ( curq == prevq ) return prevq - TENLOGTWO;*/
}

int BayesQualityGenome::AddTwoQualityValues(int curq, int prevq) {
	// TODO: rewrite correctly
	return min(curq, prevq);
	
	// there are a few shortcuts here
	/*if ( curq < prevq - 1 ) return curq;
	else if ( curq > prevq + 1 ) return prevq;
	else if ( curq == prevq + 1 ) return prevq - TENLOGELEVENTENTH;
	else if ( curq == prevq - 1 ) return curq - TENLOGELEVENTENTH;
	else if ( curq == prevq ) return prevq - TENLOGTWO;*/
}

float BayesQualityGenome::AddThreeQualityValues(float diagq, float leftq, float rightq) {
	return AddTwoQualityValues(diagq, AddTwoQualityValues(leftq, rightq));
}

int BayesQualityGenome::AddThreeQualityValues(int diagq, int leftq, int rightq) {
	return AddTwoQualityValues(diagq, AddTwoQualityValues(leftq, rightq));
}


}




/* older version without the indels
	// compute the sum of all quality values
	int totalQ = accumulate(r.second.begin(), r.second.end(), 0);
	fill(qv.begin(), qv.end(), totalQ);
	// for each element of the read
	for (size_t i = 0; i < r.first.size(); ++i) {
		// subtract its quality value if the corresponding elements match
		for (size_t j = 0; j < qv.size(); ++j) {
			if ( r.first[i] == genome_[i+j] ) {
				qv[j] -= r.second[i];
			}
		}
	} */
	
