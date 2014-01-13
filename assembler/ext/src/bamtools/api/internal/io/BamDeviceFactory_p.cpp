// ***************************************************************************
// BamDeviceFactory_p.cpp (c) 2011 Derek Barnett
// Marth Lab, Department of Biology, Boston College
// ---------------------------------------------------------------------------
// Last modified: 25 September 2011 (DB)
// ---------------------------------------------------------------------------
// Creates built-in concrete implementations of IBamIODevices
// ***************************************************************************

#include "api/internal/io/BamDeviceFactory_p.h"
#include "api/internal/io/BamFile_p.h"
using namespace BamTools;
using namespace BamTools::Internal;

#include <iostream>
using namespace std;

IBamIODevice* BamDeviceFactory::CreateDevice(const string& source) {
    // otherwise assume a "normal" file
    return new BamFile(source);
}
