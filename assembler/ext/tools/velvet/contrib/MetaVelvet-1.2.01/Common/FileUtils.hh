#ifndef _FILE_UTILS_HH_
#define _FILE_UTILS_HH_
#include <iostream>
#include <fstream>
#include <string>
#include <stdlib.h>

using namespace std;

class FileUtils {
public:
  static bool isValidFile(const string& filename){
    ifstream ifs( filename.c_str() );
    if( ifs.fail() ) return false;
    ifs.close();
    return true;
  }

  static void open(ifstream& ifs, const string& filename){
    try {
      ifs.open( filename.c_str() );
      if( ifs.fail() )
	throw (string)"[FileUtils] File open error: " + filename;
    } catch( string errorMessage ){
      cerr << errorMessage << endl;
      exit(-1);
    }
  }

  static void open(ofstream& ofs, const string& filename){
    try {
      ofs.open( filename.c_str() );
      if( ofs.fail() )
	throw (string)"[FileUtils] File open error: " + filename;
    } catch( string errorMessage ){
      cerr << errorMessage << endl;
      exit(-1);
    }
  }
};

#endif // FILE_UTILS
