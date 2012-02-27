/***************************************************************************
 * Title:          RuleList.h 
 * Author:         Mark Chaisson
 * Created:        2008
 * Last modified:  10/31/2009
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#ifndef RULE_LIST_H_
#define RULE_LIST_H_

#include <iostream>
#include <string>
#include <vector>
#include <regex.h>
using namespace std;
class Rule {
public:
	std::string regexpStr;
	ssize_t cloneLength;
	ssize_t cloneVar;
	ssize_t type;
	string forward, reverse;
	regex_t compRegex;
};


typedef std::vector<Rule> RuleList;

void ParseRuleFile(std::string &ruleFile, std::vector<Rule> &rules,
									 std::ostream &report = std::cout);

ssize_t GetReadRule(RuleList &rules, std::string &readTitle, ssize_t &readType);

#endif
