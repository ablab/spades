// Copyright (c) 2000-2003 Whitehead Institute for Biomedical Research


// Fastb: convert a fasta file into a fastb file. If FASTAMB=True, it also
// generate the relative .fastamb file.

#include "Basevector.h"
#include "Bitvector.h"
#include "MainTools.h"
#include "FetchReads.h"
#include "FetchReadsAmb.h"

int main( int argc, char *argv[] )
{
     RunTime( );

     BeginCommandArguments;
     CommandArgument_String_Doc(PRE,
			    "Folder, in which specified fasta FILE resides and into"
			    "which resulting fastb will be written. If not specified,"
                            "the value is loaded automatically from ARACHNE_PRE"
			    "environment variable");
     CommandArgument_String_Abbr_Doc(FILE, IN,
			    "Name, relative to PRE folder, of the fasta file to "
			    "convert. Must end with .fasta");
     CommandArgument_Bool_OrDefault(FASTAMB, False);
     EndCommandArguments;

     String filename;
     if ( PRE.empty() ) {
       filename = FILE;
     } else {
       if ( PRE[PRE.size()-1]!='/' ) {
	 filename = PRE+"/"+FILE;
       } else {
	 filename = PRE + FILE;
       }
     }
     ForceAssert( filename.Contains( ".fasta", -1 ) );
     String new_filename = filename.substr( 0, filename.size( ) - 6 ) + ".fastb";

     vecbasevector EE;
     FetchReads( EE, 0, filename );
     EE.WriteAll( new_filename );
     
     if ( FASTAMB ) {
       vecbitvector EEamb;
       FetchReadsAmb( EEamb, filename );
       String new_filename_amb = filename.substr ( 0, filename.size( ) - 6 ) + ".fastamb";
       EEamb.WriteAll( new_filename_amb );    
     }    
     return 0;
}
