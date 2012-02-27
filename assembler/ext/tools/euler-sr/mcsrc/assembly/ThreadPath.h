/***************************************************************************
 * Title:          ThreadPath.h 
 * Author:         Mark Chaisson
 * Created:        2008
 * Last modified:  11/26/2009
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#ifndef THREAD_PATH_H_
#define THREAD_PATH_H_

class ThreadPathInterval {
public:
	ssize_t edge;
	ssize_t length;
	ssize_t pos;
	ThreadPathInterval(ssize_t e, ssize_t l, ssize_t p) : edge(e), length(l), pos(p) {};
	ThreadPathInterval & operator=(const ThreadPathInterval &rhs) {
		if (this != &rhs) {
			edge = rhs.edge;
			length = rhs.length;
			pos = rhs.pos;
		}
		return *this;
	}
};

typedef std::list<ThreadPathInterval> ThreadPath;
typedef std::vector<ThreadPathInterval> ThreadPathVector;

#endif
