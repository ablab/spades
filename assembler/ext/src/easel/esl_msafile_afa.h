/* i/o of multiple sequence alignment files in aligned FASTA format
 */
#ifndef eslMSAFILE_AFA_INCLUDED
#define eslMSAFILE_AFA_INCLUDED
#include <esl_config.h>

#include "esl_msa.h"
#include "esl_msafile.h"

extern int esl_msafile_afa_SetInmap     (ESL_MSAFILE *afp);
extern int esl_msafile_afa_GuessAlphabet(ESL_MSAFILE *afp, int *ret_type);
extern int esl_msafile_afa_Read         (ESL_MSAFILE *afp, ESL_MSA **ret_msa);
extern int esl_msafile_afa_Write        (FILE *fp, const ESL_MSA *msa);

#endif /* eslMSAFILE_AFA_INCLUDED */

