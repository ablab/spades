#include "stdafx.h"
/*
  This file is a part of KMC software distributed under GNU GPL 3 licence.
  The homepage of the KMC project is http://sun.aei.polsl.pl/kmc
  
  Authors: Sebastian Deorowicz, Agnieszka Debudaj-Grabysz, Marek Kokot
  
  Version: 2.2.0
  Date   : 2015-04-15
*/

#include "mmer.h"


uint32 CMmer::norm5[];
uint32 CMmer::norm6[];
uint32 CMmer::norm7[];
uint32 CMmer::norm8[];

CMmer::_si CMmer::_init;


//--------------------------------------------------------------------------
CMmer::CMmer(uint32 _len)
{
	switch (_len)
	{
	case 5:
		norm = norm5;
		break;
	case 6:
		norm = norm6;
		break;
	case 7:
		norm = norm7;
		break;
	case 8:
		norm = norm8;
		break;
	default:
		break;
	}
	len = _len;
	mask = (1 << _len * 2) - 1;
	str = 0;
}

//--------------------------------------------------------------------------

