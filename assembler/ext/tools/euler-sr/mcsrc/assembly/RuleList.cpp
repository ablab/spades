/***************************************************************************
 * Title:          RuleList.cpp 
 * Author:         Mark Chaisson
 * Created:        2008
 * Last modified:  11/24/2009
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/

#include <fstream>
#include <iostream>

#include "RuleList.h"
#include "utils.h"
#include "ParseTitle.h"

using namespace std;


ssize_t IsComment(string &line) {
	ssize_t i = 0;
	while (i < line.size() and isspace(line[i])) { i++; }
	if (i < line.size() and line[i] == '#')
		return 1;
	else 
		return 0;
}
		

void ParseRuleFile(std::string &ruleFile, std::vector<Rule> &rules,
									 std::ostream &report) {
	std::ifstream ruleIn;
	openck(ruleFile,ruleIn, std::ios::in, report);
	std::string ruleLine;
	//UNUSED// ssize_t length;
	//UNUSED// ssize_t cloneVar;
	Rule rule;
	while(std::getline(ruleIn, ruleLine)) {
		if (IsComment(ruleLine)) {
			continue;
		}
		rule.regexpStr = "";
		ExtractQuotedString(ruleLine, rule.regexpStr);
		if (rule.regexpStr == "")
			break;
		ParseKeyword(ruleLine, "CloneLength", rule.cloneLength);
		ParseKeyword(ruleLine, "CloneVar", rule.cloneVar);
		if (ParseKeyword(ruleLine, "Type", rule.type) == 0) {
			rule.type = rules.size();
		}
		/*
			ParseKeyword(ruleLine, "Forward", rule.forward);
			ParseKeyword(ruleLine, "Reverse", rule.reverse);
		*/
		cout << "got rule: " << rule.regexpStr << endl;
		rules.push_back(rule);
		if (regcomp(&rules[rules.size()-1].compRegex, rule.regexpStr.c_str(), REG_EXTENDED)) {
			cout << "ERROR Parsing " << rule.regexpStr << endl;
			exit(1);
		}

	}
	ssize_t i;
	for (i= 0; i < rules.size(); i++) { 
		//		rules[i].regexp = rules[i].regexpStr;
	}
}


ssize_t GetReadRule(RuleList &rules, std::string &readTitle, ssize_t &readType) {

	ssize_t r;
	const char* readTitleStr = readTitle.c_str();
	_SZT_ nmatch = 4;
	regmatch_t pmatch[4];
	for (r = 0; r < rules.size(); r++ ){
		ssize_t mv;
		if (!(mv = regexec(&rules[r].compRegex, readTitleStr, nmatch, pmatch, 0))) {
			readType = r;
			return 1;
		}
	}
	return 0;
}
