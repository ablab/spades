/*
Copyright 2007, 2008 Daniel Zerbino (zerbino@ebi.ac.uk)

    This file is part of Velvet.

    Velvet is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    Velvet is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with Velvet; if not, write to the Free Software
    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA

*/
/****************************************************************\
*                                                                *
*  Efficient Memory Allocation Routines                          *
*                                                                *
*  Guy St.C. Slater..   mailto:guy@ebi.ac.uk                     *
*  Copyright (C) 2000-2005.  All Rights Reserved.                *
*                                                                *
*  This source code is distributed under the terms of the        *
*  GNU Lesser General Public License. See the file COPYING       *
*  or http://www.fsf.org/copyleft/lesser.html for details        *
*                                                                *
*  If you use this code, please keep this notice intact.         *
*                                                                *
\****************************************************************/

#ifndef INCLUDED_RECYCLEBIN_H
#define INCLUDED_RECYCLEBIN_H

typedef struct recycleBin_st RecycleBin;

// Constructor, Destructor
RecycleBin *newRecycleBin(size_t node_size, int nodes_per_chunk);
void destroyRecycleBin(RecycleBin * recycle_bin);

// Use
void *allocatePointer(RecycleBin * recycle_bin);
void deallocatePointer(RecycleBin * recycle_bin, void *data);

// Stats
size_t RecycleBin_memory_usage(RecycleBin * recycle_bin);
size_t recycleBinFreeSpace(RecycleBin * recycle_bin);
size_t recycleBinAvailablePointers(RecycleBin * recycle_bin);

#endif				/* INCLUDED_RECYCLEBIN_H */
