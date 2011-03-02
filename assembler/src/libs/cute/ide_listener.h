/*********************************************************************************
 * This file is part of CUTE.
 *
 * CUTE is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * CUTE is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with CUTE.  If not, see <http://www.gnu.org/licenses/>.
 *
 * Copyright 2007 Peter Sommerlad
 *
 *********************************************************************************/

#ifndef IDE_LISTENER_H_
#define IDE_LISTENER_H_
#include "eclipse_listener.h"
#include "vstudio_listener.h"
namespace cute
{
// assume that we compile with gnu when using Eclipse CDT.
// vstudio_listener is broken for VS later than 2003, TODO!
#if defined(__GNUG__)
typedef eclipse_listener ide_listener;
#else
typedef vstudio_listener ide_listener;
#endif
}

#endif /*IDE_LISTENER_H_*/
