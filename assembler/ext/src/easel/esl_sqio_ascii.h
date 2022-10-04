/* Unaligned ascii sequence file i/o.
 */
#ifndef eslSQIO_ASCII_INCLUDED
#define eslSQIO_ASCII_INCLUDED
#include "esl_config.h"

#include <stdio.h>
#ifdef HAVE_SYS_TYPES_H
#include <sys/types.h>
#endif

#include "esl_msa.h"
#include "esl_msafile.h"
#include "esl_sq.h"
#include "esl_sqio.h"

/* set the max residue count to 1 meg when reading a block */
#define MAX_RESIDUE_COUNT (1024 * 1024)

/* forward declaration */
struct esl_sqio_s;

/* ESL_SQASCII:
 * An open sequence file for reading.
 */
typedef struct esl_sqascii_s {
  
  FILE *fp;           	      /* Open file ptr                            */
  char  errbuf[eslERRBUFSIZE];/* parse error mesg.  Size must match msa.h */

  int   do_gzip;	      /* TRUE if we're reading from gzip -dc pipe */
  int   do_stdin;	      /* TRUE if we're reading from stdin         */
  int   do_buffer;            /* TRUE if we're reading from a buffer      */

  /* all input first gets buffered in memory; this gives us enough
   * recall to use Guess*() functions even in nonrewindable streams
   */
  char    *mem;		      /* buffered input                           */
  int      allocm;	      /* <mem> size, multiples of eslREADBUFSIZE  */
  int      mn;		      /* number of chars in <mem> (up to allocm)  */
  int      mpos;	      /* pos of next <buf> to load from <mem>     */
  off_t    moff;	      /* disk offset to start of <mem>            */
  int      is_recording;      /* TRUE if we need to keep buffering more   */

  /* input is either character-based [fread()] or line-based (esl_fgets())*/
  char    *buf;		      /* buffer for fread() or fgets() input      */
  off_t    boff;	      /* disk offset to start of buffer           */
  int      balloc;	      /* allocated size of buf                    */
  int      nc;		      /* #chars in buf (usually full, less at EOF)*/ 
  int      bpos;	      /* current position in the buffer (0..nc-1) */
  int64_t  L;		      /* #residues seen so far in current seq     */
  int64_t  linenumber;	      /* What line of the file  (1..N; -1=unknown)*/
  off_t    bookmark_offset;   /* bookmark fwd position before reversing...*/
  int64_t  bookmark_linenum;  /* in both linenumber and disk offset       */

  /* Format-specific configuration                                           */
  int   is_linebased;	      /* TRUE for fgets() parsers; FALSE for fread() */
  int   eof_is_ok;	      /* TRUE if record can end on EOF               */
  int  (*parse_header)(struct esl_sqio_s *, ESL_SQ *sq);
  int  (*skip_header) (struct esl_sqio_s *, ESL_SQ *sq);
  int  (*parse_end)   (struct esl_sqio_s *, ESL_SQ *sq); 

  /* MSA files can be read as sequential seq files.                     */
  ESL_MSAFILE  *afp;	      /* open ESL_MSAFILE for reading           */
  ESL_MSA      *msa;	      /* preloaded alignment to draw seqs from  */
  int           idx;	      /* index of next seq to return, 0..nseq-1 */

  /* SSI indexes allow fast random access of records in a seq file          */
  char    *ssifile;	      /* path to expected SSI index file            */
  int      rpl;		      /* residues per line in file; -1=unset 0=inval*/
  int      bpl;		      /* bytes per line in file; -1=unset, 0=inval  */
  int      currpl;	      /* residues on current line (-1=unknown)      */
  int      curbpl;	      /* bytes on current line    (-1=unknown)      */
  int      prvrpl;	      /* residues on previous line                  */
  int      prvbpl;	      /* bytes on previous line                     */
  ESL_SSI *ssi;		      /* open ESL_SSI index, or NULL if none        */
} ESL_SQASCII_DATA;


extern int  esl_sqascii_Open(char *seqfile, int format, struct esl_sqio_s *sqfp);
extern int  esl_sqascii_WriteFasta(FILE *fp, ESL_SQ *s, int update);
extern int  esl_sqascii_Parse(char *buf, int size, ESL_SQ *s, int format);

#endif /*eslSQIO_ASCII_INCLUDED*/
