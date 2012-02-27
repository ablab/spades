/***************************************************************************
 * Title:          nathans_time_hacks.h 
 * Author:         Mark Chaisson
 * Created:        2007
 * Last modified:  12/18/2008
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
/*
** Copyright 2004 Ross Lippert
** 
** This file is part of bbbwt.
** 
** bbbwt is free software; you can redistribute it and/or modify
** it under the terms of the GNU General Public License as published by
** the Free Software Foundation; either version 2 of the License, or
** (at your option) any later version.
** 
** bbbwt is distributed in the hope that it will be useful,
** but WITHOUT ANY WARRANTY; without even the implied warranty of
** MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
** GNU General Public License for more details.
** 
** You should have received a copy of the GNU General Public License
** along with bbbwt; if not, write to the Free Software
** Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
**
*/
#ifndef NATHANS_TIME_HACKS_H
#define NATHANS_TIME_HACKS_H

// Nathan's time counter
#include <time.h>
static time_t tictoctime_;
#define TIC {tictoctime_ = time(NULL);}
#define TOC(F) { time_t t; time(&t); t-=tictoctime_; fprintf(F,"Elapsed time: %d:%02d:%02d",int(t/3600),int((t%3600)/60),int(t%60)); fflush(F);} 

#endif
