/* esl_subcmd : utilities for commandline programs that take subcommands
 * 
 * See also:  
 *    esl_getopts : command line argument parsing
 */
#ifndef eslSUBCMD_INCLUDED
#define eslSUBCMD_INCLUDED
#include <esl_config.h>

#include "esl_getopts.h"

typedef struct esl_subcmd_s {
  int  (*func)(const char *topcmd, const struct esl_subcmd_s *sub, int argc, char **argv);
  char *subcmd;
  int   nargs;
  char *usage;
  char *description;
} ESL_SUBCMD;

extern ESL_GETOPTS *esl_subcmd_CreateDefaultApp(const char *topcmd, const ESL_SUBCMD *sub, const ESL_OPTIONS *suboptions, int argc, char **argv);


#endif /*eslSUBCMD_INCLUDED*/
