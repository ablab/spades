#ifndef _STDIO_HACK_H_
#define _STDIO_HACK_H_

#include <stdio.h>
#include "io_lib/mFILE.h"

#ifdef __cplusplus
extern "C" {
#endif

/*
 * This file will define many of the stdio functions to use the
 * in-memory versions here. These are sufficient to allow the abi,
 * scf, etc reading code (but not writing) to transparently operate by
 * loading the entire file into memory and then doing in-memory
 * manipulation instead of on-disk manipulation.
 *
 * The key use for this though is to allow generation of fake FILE pointers
 * so that reading from tar files does not requiring writing to disk and
 * reading back again.
 */

#define FILE mFILE
#define fopen mfopen
#define fclose mfclose
#define fseek mfseek
#define ftell mftell
#define rewind mrewind
#ifdef feof
#    undef feof
#endif
#define feof mfeof
#define fread mfread
#define fwrite mfwrite
#define fgetc mfgetc
#define ungetc mungetc
#define fgets mfgets
#define fflush mfflush
#define fprintf mfprintf

#ifdef stdin
#    undef stdin
#    undef stdout
#    undef stderr
#endif
#define stdin mstdin()
#define stdout mstdout()
#define stderr mstderr()

#define fread_abi  mfread_abi
#define fwrite_abi mfwrite_abi
#define fread_alf  mfread_alf
#define fwrite_alf mfwrite_alf
#define fread_ctf  mfread_ctf
#define fwrite_ctf mfwrite_ctf
#define fread_pln  mfread_pln
#define fwrite_pln mfwrite_pln
#define fread_scf  mfread_scf
#define fwrite_scf mfwrite_scf
#define fread_ztr  mfread_ztr
#define fwrite_ztr mfwrite_ztr

#define exp_fread_info exp_mfread_info
#define exp_print_file exp_print_mfile

#define open_trace_file open_trace_mfile
#define fread_reading mfread_reading
#define fwrite_reading mfwrite_reading

#ifdef __cplusplus
}
#endif

#endif /* _STDIO_HACK_H_ */
