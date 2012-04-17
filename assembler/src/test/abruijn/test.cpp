//***************************************************************************
//* Copyright (c) 2011-2012 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

#include "cute.h"
#include "ide_listener.h"
#include "cute_runner.h"
#include "hash_test.hpp"
#include "spellgenome_test.hpp"
#include "common/logging.hpp"

DECL_PROJECT_LOGGER("at")

void runSuite() {
	cute::suite s;
	s += HashSuite();
	s += SpellingGenomeSuite();
	cute::ide_listener lis;
	cute::makeRunner(lis)(s, "The Suite");
}

int main() {
	runSuite();
	return 0;
}
